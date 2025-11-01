"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
from __future__ import annotations
import sys
from typing import Any, Literal

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

from collections.abc import Sequence

from pydantic import Field, field_validator, model_validator

from .base import _ClassWithInit

# Type aliases - using forward references for types defined later in this file
# Note: Grid types are defined later in this file, so we use string forward references
PICMI_Grid = (
    "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | "
    "PICMI_Cartesian3DGrid | PICMI_CylindricalGrid"
)

# ---------------
# Physics objects
# ---------------

class PICMI_ElectromagneticSolver(_ClassWithInit):
    """
    Electromagnetic field solver.
    
    The advance method used to solve Maxwell's equations. The default method is code dependent.
    
    Method options:
    
    - 'Yee': standard solver using the staggered Yee grid (https://doi.org/10.1109/TAP.1966.1138693)
    - 'CKC': solver with the extended Cole-Karkkainen-Cowan stencil with better dispersion properties (https://doi.org/10.1103/PhysRevSTAB.16.041303)
    - 'Lehe': CKC-style solver with modified dispersion (https://doi.org/10.1103/PhysRevSTAB.16.021301)
    - 'PSTD': Spectral solver with finite difference in time domain, e.g., Q. H. Liu, Letters 15 (3) (1997) 158–165
    - 'PSATD': Spectral solver with analytic in time domain (https://doi.org/10.1016/j.jcp.2013.03.010)
    - 'DS': Directional Splitting after Yasuhiko Sentoku (https://doi.org/10.1140/epjd/e2014-50162-y)
    - 'ECT': Enlarged Cell Technique solver, allowing internal conductors (https://doi.org/10.1109/APS.2005.1551259)
    """
    methods_list: list[str] = ['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', 'ECT']
    
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid object for the diagnostic")
    method: Literal['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', 'ECT'] | None = Field(
        default=None,
        description="The advance method use to solve Maxwell's equations. The default method is code dependent."
    )
    stencil_order: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers. Order of stencil for each axis (-1=infinite)"
    )
    cfl: float | None = Field(
        default=None,
        description="Fraction of the Courant-Friedrich-Lewy criteria [1]"
    )
    source_smoother: "PICMI_BinomialSmoother | None" = Field(
        default=None,
        description="Smoother instance. Smoother object to apply to the sources"
    )
    field_smoother: "PICMI_BinomialSmoother | None" = Field(
        default=None,
        description="Smoother instance. Smoother object to apply to the fields"
    )
    subcycling: int | None = Field(
        default=None,
        description="Level of subcycling for the GPSTD solver"
    )
    galilean_velocity: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats. Velocity of Galilean reference frame [m/s]"
    )
    divE_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(E) cleaning if True"
    )
    divB_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(B) cleaning if True"
    )
    pml_divE_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(E) cleaning in the PML if True"
    )
    pml_divB_cleaning: bool | None = Field(
        default=None,
        description="Solver uses div(B) cleaning in the PML if True"
    )
    
    @field_validator('method')
    @classmethod
    def _validate_method(cls, v):
        if v is not None and v not in PICMI_ElectromagneticSolver.methods_list:
            raise ValueError(f'method must be one of {", ".join(PICMI_ElectromagneticSolver.methods_list)}')
        return v


class PICMI_ElectrostaticSolver(_ClassWithInit):
    """
    Electrostatic field solver.
    """
    methods_list: list[str] = ['FFT', 'Multigrid']
    
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid instance. Grid object for the diagnostic")
    method: Literal['FFT', 'Multigrid'] | None = Field(
        default=None,
        description="String. One of 'FFT', or 'Multigrid'"
    )
    required_precision: float | None = Field(
        default=None,
        description="Level of precision required for iterative solvers"
    )
    maximum_iterations: int | None = Field(
        default=None,
        description="Maximum number of iterations for iterative solvers"
    )
    
    @field_validator('method')
    @classmethod
    def _validate_method(cls, v):
        if v is not None and v not in PICMI_ElectrostaticSolver.methods_list:
            raise ValueError(f'method must be one of {", ".join(PICMI_ElectrostaticSolver.methods_list)}')
        return v


class PICMI_MagnetostaticSolver(_ClassWithInit):
    """
    Magnetostatic field solver.
    """
    methods_list: list[str] = ['FFT', 'Multigrid']
    
    grid: "PICMI_Cartesian1DGrid | PICMI_Cartesian2DGrid | PICMI_Cartesian3DGrid | PICMI_CylindricalGrid" = Field(description="Grid instance. Grid object for the diagnostic")
    method: Literal['FFT', 'Multigrid'] | None = Field(
        default=None,
        description="String. One of 'FFT', or 'Multigrid'"
    )
    
    @field_validator('method')
    @classmethod
    def _validate_method(cls, v):
        if v is not None and v not in PICMI_MagnetostaticSolver.methods_list:
            raise ValueError(f'method must be one of {", ".join(PICMI_MagnetostaticSolver.methods_list)}')
        return v


# ------------------
# Numeric Objects
# ------------------


class PICMI_BinomialSmoother(_ClassWithInit):
    """
    Describes a binomial smoother operator (applied to grids).
    """
    n_pass: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers. Number of passes along each axis"
    )
    compensation: Sequence[bool] | None = Field(
        default=None,
        description="Vector of booleans, optional. Flags whether to apply compensation along each axis"
    )
    stride: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Stride along each axis"
    )
    alpha: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats, optional. Smoothing coefficients along each axis"
    )


class PICMI_Cartesian1DGrid(_ClassWithInit):
    """
    One-dimensional Cartesian grid.
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    number_of_dimensions: int = 1

    # Vector form parameters (preferred)
    number_of_cells: Sequence[int] | None = Field(
        default=None,
        min_length=1,
        max_length=1,
        description="Vector of integers. Number of cells along each axis (number of nodes is number_of_cells+1). Either this or nx must be specified (not both)."
    )
    lower_bound: Sequence[float] | None = Field(
        default=None,
        min_length=1,
        max_length=1,
        description="Vector of floats. Position of the node at the lower bound [m]. Either this or xmin must be specified (not both)."
    )
    upper_bound: Sequence[float] | None = Field(
        default=None,
        min_length=1,
        max_length=1,
        description="Vector of floats. Position of the node at the upper bound [m]. Either this or xmax must be specified (not both)."
    )
    lower_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=1,
        max_length=1,
        description="Vector of strings. Conditions at lower boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_xmin must be specified (not both)."
    )
    upper_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=1,
        max_length=1,
        description="Vector of strings. Conditions at upper boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_xmax must be specified (not both)."
    )

    # Individual named parameters (alternative form)
    nx: int | None = Field(
        default=None,
        description="Integer. Number of cells along X (number of nodes=nx+1). Either this or number_of_cells must be specified (not both)."
    )
    xmin: float | None = Field(
        default=None,
        description="Float. Position of first node along X [m]. Either this or lower_bound must be specified (not both)."
    )
    xmax: float | None = Field(
        default=None,
        description="Float. Position of last node along X [m]. Either this or upper_bound must be specified (not both)."
    )
    bc_xmin: str | None = Field(
        default=None,
        description="String. Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or lower_boundary_conditions must be specified (not both)."
    )
    bc_xmax: str | None = Field(
        default=None,
        description="String. Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or upper_boundary_conditions must be specified (not both)."
    )

    # Particle boundary parameters (vector form)
    lower_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle lower bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    upper_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle upper bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    lower_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at lower boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )
    upper_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at upper boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )

    # Particle boundary parameters (individual form)
    xmin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along X [m]."
    )
    xmax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along X [m]."
    )
    bc_xmin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_xmax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal."
    )

    # Other parameters
    moving_window_velocity: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats, optional. Moving frame velocity [m/s]"
    )
    refined_regions: list[list] = Field(
        default_factory=list,
        description="List of lists, optional. List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level (1 being first, 2 being second etc), lo and hi being vectors of length 1 specifying the extent of the region, and refinement_factor defaulting to [2] (relative to next lower level)."
    )
    guard_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of guard cells used along each direction"
    )
    pml_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of Perfectly Matched Layer (PML) cells along each direction"
    )

    @model_validator(mode='before')
    @classmethod
    def _normalize_parameters(cls, data: Any) -> Any:
        """Normalize individual parameters to vector form before validation"""
        if not isinstance(data, dict):
            return data
        
        data = dict(data)
        
        # Validate either/or constraints for grid parameters
        number_of_cells = data.get('number_of_cells')
        nx = data.get('nx')
        if (number_of_cells is None) == (nx is None):
            raise ValueError('Exactly one of number_of_cells or nx must be specified')
        
        lower_bound = data.get('lower_bound')
        xmin = data.get('xmin')
        if (lower_bound is None) == (xmin is None):
            raise ValueError('Exactly one of lower_bound or xmin must be specified')
        
        upper_bound = data.get('upper_bound')
        xmax = data.get('xmax')
        if (upper_bound is None) == (xmax is None):
            raise ValueError('Exactly one of upper_bound or xmax must be specified')
        
        lower_boundary_conditions = data.get('lower_boundary_conditions')
        bc_xmin = data.get('bc_xmin')
        if (lower_boundary_conditions is None) == (bc_xmin is None):
            raise ValueError('Exactly one of lower_boundary_conditions or bc_xmin must be specified')
        
        upper_boundary_conditions = data.get('upper_boundary_conditions')
        bc_xmax = data.get('bc_xmax')
        if (upper_boundary_conditions is None) == (bc_xmax is None):
            raise ValueError('Exactly one of upper_boundary_conditions or bc_xmax must be specified')
        
        # Convert individual parameters to vectors if vectors not provided
        # If both provided, vector takes precedence (as per original docstring)
        if number_of_cells is None and nx is not None:
            data['number_of_cells'] = [nx]
            data.pop('nx', None)
        elif number_of_cells is not None and nx is not None:
            data.pop('nx', None)  # Prefer vector form
        
        if lower_bound is None and xmin is not None:
            data['lower_bound'] = [xmin]
            data.pop('xmin', None)
        elif lower_bound is not None and xmin is not None:
            data.pop('xmin', None)
        
        if upper_bound is None and xmax is not None:
            data['upper_bound'] = [xmax]
            data.pop('xmax', None)
        elif upper_bound is not None and xmax is not None:
            data.pop('xmax', None)
        
        if lower_boundary_conditions is None and bc_xmin is not None:
            data['lower_boundary_conditions'] = [bc_xmin]
            data.pop('bc_xmin', None)
        elif lower_boundary_conditions is not None and bc_xmin is not None:
            data.pop('bc_xmin', None)
        
        if upper_boundary_conditions is None and bc_xmax is not None:
            data['upper_boundary_conditions'] = [bc_xmax]
            data.pop('bc_xmax', None)
        elif upper_boundary_conditions is not None and bc_xmax is not None:
            data.pop('bc_xmax', None)
        
        # Handle particle boundaries
        if data.get('lower_bound_particles') is None:
            if data.get('xmin_particles') is not None:
                data['lower_bound_particles'] = [data.pop('xmin_particles')]
        if data.get('upper_bound_particles') is None:
            if data.get('xmax_particles') is not None:
                data['upper_bound_particles'] = [data.pop('xmax_particles')]
        if data.get('lower_boundary_conditions_particles') is None:
            if data.get('bc_xmin_particles') is not None:
                data['lower_boundary_conditions_particles'] = [data.pop('bc_xmin_particles')]
        if data.get('upper_boundary_conditions_particles') is None:
            if data.get('bc_xmax_particles') is not None:
                data['upper_boundary_conditions_particles'] = [data.pop('bc_xmax_particles')]
        
        return data
    
    @model_validator(mode='after')
    def _validate_and_normalize(self) -> Self:
        """Validate dimensions and normalize particle boundaries"""
        
        # Handle particle boundaries - default to field boundaries if not specified
        if self.lower_bound_particles is None:
            if self.xmin_particles is None:
                self.lower_bound_particles = list(self.lower_bound)
            else:
                self.lower_bound_particles = [self.xmin_particles]
                self.xmin_particles = None
        if self.upper_bound_particles is None:
            if self.xmax_particles is None:
                self.upper_bound_particles = list(self.upper_bound)
            else:
                self.upper_bound_particles = [self.xmax_particles]
                self.xmax_particles = None
        if self.lower_boundary_conditions_particles is None:
            if self.bc_xmin_particles is None:
                self.lower_boundary_conditions_particles = list(self.lower_boundary_conditions)
            else:
                self.lower_boundary_conditions_particles = [self.bc_xmin_particles]
                self.bc_xmin_particles = None
        if self.upper_boundary_conditions_particles is None:
            if self.bc_xmax_particles is None:
                self.upper_boundary_conditions_particles = list(self.upper_boundary_conditions)
            else:
                self.upper_boundary_conditions_particles = [self.bc_xmax_particles]
                self.bc_xmax_particles = None
        
        # Dimensions are validated by Field(min_length=1, max_length=1) constraints
        
        # Process refined regions
        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2])
            if len(region[1]) != 1:
                raise ValueError('The lo extent of the refined region must be a vector of length 1')
            if len(region[2]) != 1:
                raise ValueError('The hi extent of the refined region must be a vector of length 1')
            if len(region[3]) != 1:
                raise ValueError('The refinement factor of the refined region must be a vector of length 1')
        
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_CylindricalGrid(_ClassWithInit):
    """
    Axisymmetric, cylindrical grid.
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    number_of_dimensions: int = 2

    # Vector form parameters (preferred)
    number_of_cells: Sequence[int] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of integers. Number of cells along each axis (number of nodes is number_of_cells+1). Either this or nr and nz must be specified (not both)."
    )
    lower_bound: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats. Position of the node at the lower bound [m]. Either this or rmin and zmin must be specified (not both)."
    )
    upper_bound: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats. Position of the node at the upper bound [m]. Either this or rmax and zmax must be specified (not both)."
    )
    lower_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings. Conditions at lower boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_rmin and bc_zmin must be specified (not both). Note: bc_rmin may be None since it will usually be the axis."
    )
    upper_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings. Conditions at upper boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_rmax and bc_zmax must be specified (not both)."
    )

    # Individual named parameters (alternative form)
    nr: int | None = Field(
        default=None,
        description="Integer. Number of cells along R (number of nodes=nr+1). Either this (along with nz) or number_of_cells must be specified (not both)."
    )
    nz: int | None = Field(
        default=None,
        description="Integer. Number of cells along Z (number of nodes=nz+1). Either this (along with nr) or number_of_cells must be specified (not both)."
    )
    n_azimuthal_modes: int | None = Field(
        default=None,
        description="Integer. Number of azimuthal modes"
    )
    rmin: float | None = Field(
        default=None,
        description="Float. Position of first node along R [m]. Either this (along with zmin) or lower_bound must be specified (not both)."
    )
    rmax: float | None = Field(
        default=None,
        description="Float. Position of last node along R [m]. Either this (along with zmax) or upper_bound must be specified (not both)."
    )
    zmin: float | None = Field(
        default=None,
        description="Float. Position of first node along Z [m]. Either this (along with rmin) or lower_bound must be specified (not both)."
    )
    zmax: float | None = Field(
        default=None,
        description="Float. Position of last node along Z [m]. Either this (along with rmax) or upper_bound must be specified (not both)."
    )
    bc_rmin: str | None = Field(
        default=None,
        description="String. Boundary condition at min R: One of open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_zmin) or lower_boundary_conditions must be specified (not both). May be None since it will usually be the axis."
    )
    bc_rmax: str | None = Field(
        default=None,
        description="String. Boundary condition at max R: One of open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_zmax) or upper_boundary_conditions must be specified (not both)."
    )
    bc_zmin: str | None = Field(
        default=None,
        description="String. Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_rmin) or lower_boundary_conditions must be specified (not both)."
    )
    bc_zmax: str | None = Field(
        default=None,
        description="String. Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_rmax) or upper_boundary_conditions must be specified (not both)."
    )

    # Particle boundary parameters (vector form)
    lower_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle lower bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    upper_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle upper bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    lower_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at lower boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )
    upper_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at upper boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )

    # Particle boundary parameters (individual form)
    rmin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along R [m]."
    )
    rmax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along R [m]."
    )
    zmin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along Z [m]."
    )
    zmax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along Z [m]."
    )
    bc_rmin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min R for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_rmax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max R for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_zmin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_zmax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal."
    )

    # Other parameters
    moving_window_velocity: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats, optional. Moving frame velocity [m/s]"
    )
    refined_regions: list[list] = Field(
        default_factory=list,
        description="List of lists, optional. List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level (1 being first, 2 being second etc), lo and hi being vectors of length 2 specifying the extent of the region, and refinement_factor defaulting to [2,2] (relative to next lower level)."
    )
    guard_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of guard cells used along each direction"
    )
    pml_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of Perfectly Matched Layer (PML) cells along each direction"
    )

    @model_validator(mode='before')
    @classmethod
    def _normalize_parameters(cls, data: Any) -> Any:
        """Normalize individual parameters to vector form before validation"""
        if not isinstance(data, dict):
            return data
        
        data = dict(data)
        
        # Validate either/or constraints for grid parameters
        number_of_cells = data.get('number_of_cells')
        nr = data.get('nr')
        nz = data.get('nz')
        if (number_of_cells is None) == ((nr is not None) and (nz is not None)):
            raise ValueError('Exactly one of number_of_cells or (nr and nz) must be specified')
        
        lower_bound = data.get('lower_bound')
        rmin = data.get('rmin')
        zmin = data.get('zmin')
        if (lower_bound is None) == ((rmin is not None) and (zmin is not None)):
            raise ValueError('Exactly one of lower_bound or (rmin and zmin) must be specified')
        
        upper_bound = data.get('upper_bound')
        rmax = data.get('rmax')
        zmax = data.get('zmax')
        if (upper_bound is None) == ((rmax is not None) and (zmax is not None)):
            raise ValueError('Exactly one of upper_bound or (rmax and zmax) must be specified')
        
        lower_boundary_conditions = data.get('lower_boundary_conditions')
        bc_rmin = data.get('bc_rmin')
        bc_zmin = data.get('bc_zmin')
        # Special case: bc_rmin can be None (axis), so only check bc_zmin
        if (lower_boundary_conditions is None) and (bc_zmin is None):
            raise ValueError('Either lower_boundary_conditions or bc_zmin (and optionally bc_rmin) must be specified')
        if (lower_boundary_conditions is not None) and (bc_rmin is not None or bc_zmin is not None):
            raise ValueError('Either lower_boundary_conditions or (bc_rmin and bc_zmin) must be specified (not both)')
        
        upper_boundary_conditions = data.get('upper_boundary_conditions')
        bc_rmax = data.get('bc_rmax')
        bc_zmax = data.get('bc_zmax')
        if (upper_boundary_conditions is None) == ((bc_rmax is not None) and (bc_zmax is not None)):
            raise ValueError('Exactly one of upper_boundary_conditions or (bc_rmax and bc_zmax) must be specified')
        
        # Convert individual parameters to vectors if vectors not provided
        if number_of_cells is None and nr is not None and nz is not None:
            data['number_of_cells'] = [nr, nz]
            data.pop('nr', None)
            data.pop('nz', None)
        elif number_of_cells is not None:
            data.pop('nr', None)
            data.pop('nz', None)
        
        if lower_bound is None and rmin is not None and zmin is not None:
            data['lower_bound'] = [rmin, zmin]
            data.pop('rmin', None)
            data.pop('zmin', None)
        elif lower_bound is not None:
            data.pop('rmin', None)
            data.pop('zmin', None)
        
        if upper_bound is None and rmax is not None and zmax is not None:
            data['upper_bound'] = [rmax, zmax]
            data.pop('rmax', None)
            data.pop('zmax', None)
        elif upper_bound is not None:
            data.pop('rmax', None)
            data.pop('zmax', None)
        
        if lower_boundary_conditions is None:
            bc_rmin_val = data.get('bc_rmin')
            bc_zmin_val = data.get('bc_zmin')
            if bc_zmin_val is not None:
                data['lower_boundary_conditions'] = [bc_rmin_val, bc_zmin_val]
                data.pop('bc_rmin', None)
                data.pop('bc_zmin', None)
        else:
            data.pop('bc_rmin', None)
            data.pop('bc_zmin', None)
        
        if upper_boundary_conditions is None and bc_rmax is not None and bc_zmax is not None:
            data['upper_boundary_conditions'] = [bc_rmax, bc_zmax]
            data.pop('bc_rmax', None)
            data.pop('bc_zmax', None)
        elif upper_boundary_conditions is not None:
            data.pop('bc_rmax', None)
            data.pop('bc_zmax', None)
        
        # Handle particle boundaries
        if data.get('lower_bound_particles') is None:
            rmin_p = data.get('rmin_particles')
            zmin_p = data.get('zmin_particles')
            if rmin_p is not None or zmin_p is not None:
                data['lower_bound_particles'] = [rmin_p, zmin_p]
                data.pop('rmin_particles', None)
                data.pop('zmin_particles', None)
        if data.get('upper_bound_particles') is None:
            rmax_p = data.get('rmax_particles')
            zmax_p = data.get('zmax_particles')
            if rmax_p is not None or zmax_p is not None:
                data['upper_bound_particles'] = [rmax_p, zmax_p]
                data.pop('rmax_particles', None)
                data.pop('zmax_particles', None)
        if data.get('lower_boundary_conditions_particles') is None:
            bc_rmin_p = data.get('bc_rmin_particles')
            bc_zmin_p = data.get('bc_zmin_particles')
            if bc_rmin_p is not None or bc_zmin_p is not None:
                data['lower_boundary_conditions_particles'] = [bc_rmin_p, bc_zmin_p]
                data.pop('bc_rmin_particles', None)
                data.pop('bc_zmin_particles', None)
        if data.get('upper_boundary_conditions_particles') is None:
            bc_rmax_p = data.get('bc_rmax_particles')
            bc_zmax_p = data.get('bc_zmax_particles')
            if bc_rmax_p is not None or bc_zmax_p is not None:
                data['upper_boundary_conditions_particles'] = [bc_rmax_p, bc_zmax_p]
                data.pop('bc_rmax_particles', None)
                data.pop('bc_zmax_particles', None)
        
        return data
    
    @model_validator(mode='after')
    def _validate_and_normalize(self) -> Self:
        """Validate dimensions and normalize particle boundaries"""
        # At this point, vectors should already be set by mode='before' validator
        # Handle particle boundaries - default to field boundaries if not specified
        if self.lower_bound_particles is None:
            if self.rmin_particles is None and self.zmin_particles is None:
                self.lower_bound_particles = list(self.lower_bound)
            else:
                self.lower_bound_particles = [self.rmin_particles, self.zmin_particles]
                self.rmin_particles = None
                self.zmin_particles = None
        if self.upper_bound_particles is None:
            if self.rmax_particles is None and self.zmax_particles is None:
                self.upper_bound_particles = list(self.upper_bound)
            else:
                self.upper_bound_particles = [self.rmax_particles, self.zmax_particles]
                self.rmax_particles = None
                self.zmax_particles = None
        if self.lower_boundary_conditions_particles is None:
            if self.bc_rmin_particles is None and self.bc_zmin_particles is None:
                self.lower_boundary_conditions_particles = list(self.lower_boundary_conditions)
            else:
                self.lower_boundary_conditions_particles = [self.bc_rmin_particles, self.bc_zmin_particles]
                self.bc_rmin_particles = None
                self.bc_zmin_particles = None
        if self.upper_boundary_conditions_particles is None:
            if self.bc_rmax_particles is None and self.bc_zmax_particles is None:
                self.upper_boundary_conditions_particles = list(self.upper_boundary_conditions)
            else:
                self.upper_boundary_conditions_particles = [self.bc_rmax_particles, self.bc_zmax_particles]
                self.bc_rmax_particles = None
                self.bc_zmax_particles = None
        
        # Dimensions are validated by Field(min_length=2, max_length=2) constraints
        
        # Process refined regions
        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2, 2])
            if len(region[1]) != 2:
                raise ValueError('The lo extent of the refined region must be a vector of length 2')
            if len(region[2]) != 2:
                raise ValueError('The hi extent of the refined region must be a vector of length 2')
            if len(region[3]) != 2:
                raise ValueError('The refinement factor of the refined region must be a vector of length 2')
        
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2, 2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_Cartesian2DGrid(_ClassWithInit):
    """
    Two dimensional Cartesian grid.
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    number_of_dimensions: int = 2

    # Vector form parameters (preferred)
    number_of_cells: Sequence[int] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of integers. Number of cells along each axis (number of nodes is number_of_cells+1). Either this or nx and ny must be specified (not both)."
    )
    lower_bound: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats. Position of the node at the lower bound [m]. Either this or xmin and ymin must be specified (not both)."
    )
    upper_bound: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats. Position of the node at the upper bound [m]. Either this or xmax and ymax must be specified (not both)."
    )
    lower_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings. Conditions at lower boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_xmin and bc_ymin must be specified (not both)."
    )
    upper_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings. Conditions at upper boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_xmax and bc_ymax must be specified (not both)."
    )

    # Individual named parameters (alternative form)
    nx: int | None = Field(
        default=None,
        description="Integer. Number of cells along X (number of nodes=nx+1). Either this (along with ny) or number_of_cells must be specified (not both)."
    )
    ny: int | None = Field(
        default=None,
        description="Integer. Number of cells along Y (number of nodes=ny+1). Either this (along with nx) or number_of_cells must be specified (not both)."
    )
    xmin: float | None = Field(
        default=None,
        description="Float. Position of first node along X [m]. Either this (along with ymin) or lower_bound must be specified (not both)."
    )
    xmax: float | None = Field(
        default=None,
        description="Float. Position of last node along X [m]. Either this (along with ymax) or upper_bound must be specified (not both)."
    )
    ymin: float | None = Field(
        default=None,
        description="Float. Position of first node along Y [m]. Either this (along with xmin) or lower_bound must be specified (not both)."
    )
    ymax: float | None = Field(
        default=None,
        description="Float. Position of last node along Y [m]. Either this (along with xmax) or upper_bound must be specified (not both)."
    )
    bc_xmin: str | None = Field(
        default=None,
        description="String. Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_ymin) or lower_boundary_conditions must be specified (not both)."
    )
    bc_xmax: str | None = Field(
        default=None,
        description="String. Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_ymax) or upper_boundary_conditions must be specified (not both)."
    )
    bc_ymin: str | None = Field(
        default=None,
        description="String. Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_xmin) or lower_boundary_conditions must be specified (not both)."
    )
    bc_ymax: str | None = Field(
        default=None,
        description="String. Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_xmax) or upper_boundary_conditions must be specified (not both)."
    )

    # Particle boundary parameters (vector form)
    lower_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle lower bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    upper_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle upper bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    lower_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at lower boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )
    upper_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at upper boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )

    # Particle boundary parameters (individual form)
    xmin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along X [m]."
    )
    xmax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along X [m]."
    )
    ymin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along Y [m]."
    )
    ymax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along Y [m]."
    )
    bc_xmin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_xmax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_ymin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_ymax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal."
    )

    # Other parameters
    moving_window_velocity: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats, optional. Moving frame velocity [m/s]"
    )
    refined_regions: list[list] = Field(
        default_factory=list,
        description="List of lists, optional. List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level (1 being first, 2 being second etc), lo and hi being vectors of length 2 specifying the extent of the region, and refinement_factor defaulting to [2,2] (relative to next lower level)."
    )
    guard_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of guard cells used along each direction"
    )
    pml_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of Perfectly Matched Layer (PML) cells along each direction"
    )

    @model_validator(mode='before')
    @classmethod
    def _normalize_parameters(cls, data: Any) -> Any:
        """Normalize individual parameters to vector form before validation"""
        if not isinstance(data, dict):
            return data
        
        data = dict(data)
        
        # Validate either/or constraints for grid parameters
        number_of_cells = data.get('number_of_cells')
        nx = data.get('nx')
        ny = data.get('ny')
        if (number_of_cells is None) == ((nx is not None) and (ny is not None)):
            raise ValueError('Exactly one of number_of_cells or (nx and ny) must be specified')
        
        lower_bound = data.get('lower_bound')
        xmin = data.get('xmin')
        ymin = data.get('ymin')
        if (lower_bound is None) == ((xmin is not None) and (ymin is not None)):
            raise ValueError('Exactly one of lower_bound or (xmin and ymin) must be specified')
        
        upper_bound = data.get('upper_bound')
        xmax = data.get('xmax')
        ymax = data.get('ymax')
        if (upper_bound is None) == ((xmax is not None) and (ymax is not None)):
            raise ValueError('Exactly one of upper_bound or (xmax and ymax) must be specified')
        
        lower_boundary_conditions = data.get('lower_boundary_conditions')
        bc_xmin = data.get('bc_xmin')
        bc_ymin = data.get('bc_ymin')
        if (lower_boundary_conditions is None) == ((bc_xmin is not None) and (bc_ymin is not None)):
            raise ValueError('Exactly one of lower_boundary_conditions or (bc_xmin and bc_ymin) must be specified')
        
        upper_boundary_conditions = data.get('upper_boundary_conditions')
        bc_xmax = data.get('bc_xmax')
        bc_ymax = data.get('bc_ymax')
        if (upper_boundary_conditions is None) == ((bc_xmax is not None) and (bc_ymax is not None)):
            raise ValueError('Exactly one of upper_boundary_conditions or (bc_xmax and bc_ymax) must be specified')
        
        # Convert individual parameters to vectors if vectors not provided
        if number_of_cells is None and nx is not None and ny is not None:
            data['number_of_cells'] = [nx, ny]
            data.pop('nx', None)
            data.pop('ny', None)
        elif number_of_cells is not None:
            data.pop('nx', None)
            data.pop('ny', None)
        
        if lower_bound is None and xmin is not None and ymin is not None:
            data['lower_bound'] = [xmin, ymin]
            data.pop('xmin', None)
            data.pop('ymin', None)
        elif lower_bound is not None:
            data.pop('xmin', None)
            data.pop('ymin', None)
        
        if upper_bound is None and xmax is not None and ymax is not None:
            data['upper_bound'] = [xmax, ymax]
            data.pop('xmax', None)
            data.pop('ymax', None)
        elif upper_bound is not None:
            data.pop('xmax', None)
            data.pop('ymax', None)
        
        if lower_boundary_conditions is None and bc_xmin is not None and bc_ymin is not None:
            data['lower_boundary_conditions'] = [bc_xmin, bc_ymin]
            data.pop('bc_xmin', None)
            data.pop('bc_ymin', None)
        elif lower_boundary_conditions is not None:
            data.pop('bc_xmin', None)
            data.pop('bc_ymin', None)
        
        if upper_boundary_conditions is None and bc_xmax is not None and bc_ymax is not None:
            data['upper_boundary_conditions'] = [bc_xmax, bc_ymax]
            data.pop('bc_xmax', None)
            data.pop('bc_ymax', None)
        elif upper_boundary_conditions is not None:
            data.pop('bc_xmax', None)
            data.pop('bc_ymax', None)
        
        # Handle particle boundaries
        if data.get('lower_bound_particles') is None:
            xmin_p = data.get('xmin_particles')
            ymin_p = data.get('ymin_particles')
            if xmin_p is not None or ymin_p is not None:
                data['lower_bound_particles'] = [xmin_p, ymin_p]
                data.pop('xmin_particles', None)
                data.pop('ymin_particles', None)
        if data.get('upper_bound_particles') is None:
            xmax_p = data.get('xmax_particles')
            ymax_p = data.get('ymax_particles')
            if xmax_p is not None or ymax_p is not None:
                data['upper_bound_particles'] = [xmax_p, ymax_p]
                data.pop('xmax_particles', None)
                data.pop('ymax_particles', None)
        if data.get('lower_boundary_conditions_particles') is None:
            bc_xmin_p = data.get('bc_xmin_particles')
            bc_ymin_p = data.get('bc_ymin_particles')
            if bc_xmin_p is not None or bc_ymin_p is not None:
                data['lower_boundary_conditions_particles'] = [bc_xmin_p, bc_ymin_p]
                data.pop('bc_xmin_particles', None)
                data.pop('bc_ymin_particles', None)
        if data.get('upper_boundary_conditions_particles') is None:
            bc_xmax_p = data.get('bc_xmax_particles')
            bc_ymax_p = data.get('bc_ymax_particles')
            if bc_xmax_p is not None or bc_ymax_p is not None:
                data['upper_boundary_conditions_particles'] = [bc_xmax_p, bc_ymax_p]
                data.pop('bc_xmax_particles', None)
                data.pop('bc_ymax_particles', None)
        
        return data
    
    @model_validator(mode='after')
    def _validate_and_normalize(self) -> Self:
        """Validate dimensions and normalize particle boundaries"""
        # At this point, vectors should already be set by mode='before' validator
        # Handle particle boundaries - default to field boundaries if not specified
        if self.lower_bound_particles is None:
            if self.xmin_particles is None and self.ymin_particles is None:
                self.lower_bound_particles = list(self.lower_bound)
            else:
                self.lower_bound_particles = [self.xmin_particles, self.ymin_particles]
                self.xmin_particles = None
                self.ymin_particles = None
        if self.upper_bound_particles is None:
            if self.xmax_particles is None and self.ymax_particles is None:
                self.upper_bound_particles = list(self.upper_bound)
            else:
                self.upper_bound_particles = [self.xmax_particles, self.ymax_particles]
                self.xmax_particles = None
                self.ymax_particles = None
        if self.lower_boundary_conditions_particles is None:
            if self.bc_xmin_particles is None and self.bc_ymin_particles is None:
                self.lower_boundary_conditions_particles = list(self.lower_boundary_conditions)
            else:
                self.lower_boundary_conditions_particles = [self.bc_xmin_particles, self.bc_ymin_particles]
                self.bc_xmin_particles = None
                self.bc_ymin_particles = None
        if self.upper_boundary_conditions_particles is None:
            if self.bc_xmax_particles is None and self.bc_ymax_particles is None:
                self.upper_boundary_conditions_particles = list(self.upper_boundary_conditions)
            else:
                self.upper_boundary_conditions_particles = [self.bc_xmax_particles, self.bc_ymax_particles]
                self.bc_xmax_particles = None
                self.bc_ymax_particles = None
        
        # Dimensions are validated by Field(min_length=2, max_length=2) constraints
        
        # Process refined regions
        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2, 2])
            if len(region[1]) != 2:
                raise ValueError('The lo extent of the refined region must be a vector of length 2')
            if len(region[2]) != 2:
                raise ValueError('The hi extent of the refined region must be a vector of length 2')
            if len(region[3]) != 2:
                raise ValueError('The refinement factor of the refined region must be a vector of length 2')
        
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2, 2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_Cartesian3DGrid(_ClassWithInit):
    """
    Three dimensional Cartesian grid.
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    References
    ----------
    absorbing_silver_mueller: A local absorbing boundary condition that works best under normal incidence angle.
    Based on the Silver-Mueller Radiation Condition, e.g., in

    * A. K. Belhora and L. Pichon, "Maybe Efficient Absorbing Boundary Conditions for the Finite Element Solution of 3D Scattering Problems," 1995,
      https://doi.org/10.1109/20.376322
    * B Engquist and A. Majdat, "Absorbing boundary conditions for numerical simulation of waves," 1977,
      https://doi.org/10.1073/pnas.74.5.1765
    * R. Lehe, "Electromagnetic wave propagation in Particle-In-Cell codes," 2016,
      US Particle Accelerator School (USPAS) Summer Session, Self-Consistent Simulations of Beam and Plasma Systems
      https://people.nscl.msu.edu/~lund/uspas/scs_2016/lec_adv/A1b_EM_Waves.pdf
    """
    number_of_dimensions: int = 3

    # Vector form parameters (preferred)
    number_of_cells: Sequence[int] | None = Field(
        default=None,
        min_length=3,
        max_length=3,
        description="Vector of integers. Number of cells along each axis (number of nodes is number_of_cells+1). Either this or nx, ny, and nz must be specified (not both)."
    )
    lower_bound: Sequence[float] | None = Field(
        default=None,
        min_length=3,
        max_length=3,
        description="Vector of floats. Position of the node at the lower bound [m]. Either this or xmin, ymin, and zmin must be specified (not both)."
    )
    upper_bound: Sequence[float] | None = Field(
        default=None,
        min_length=3,
        max_length=3,
        description="Vector of floats. Position of the node at the upper bound [m]. Either this or xmax, ymax, and zmax must be specified (not both)."
    )
    lower_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=3,
        max_length=3,
        description="Vector of strings. Conditions at lower boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_xmin, bc_ymin, and bc_zmin must be specified (not both)."
    )
    upper_boundary_conditions: Sequence[str] | None = Field(
        default=None,
        min_length=3,
        max_length=3,
        description="Vector of strings. Conditions at upper boundaries: periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this or bc_xmax, bc_ymax, and bc_zmax must be specified (not both)."
    )

    # Individual named parameters (alternative form)
    nx: int | None = Field(
        default=None,
        description="Integer. Number of cells along X (number of nodes=nx+1). Either this (along with ny and nz) or number_of_cells must be specified (not both)."
    )
    ny: int | None = Field(
        default=None,
        description="Integer. Number of cells along Y (number of nodes=ny+1). Either this (along with nx and nz) or number_of_cells must be specified (not both)."
    )
    nz: int | None = Field(
        default=None,
        description="Integer. Number of cells along Z (number of nodes=nz+1). Either this (along with nx and ny) or number_of_cells must be specified (not both)."
    )
    xmin: float | None = Field(
        default=None,
        description="Float. Position of first node along X [m]. Either this (along with ymin and zmin) or lower_bound must be specified (not both)."
    )
    xmax: float | None = Field(
        default=None,
        description="Float. Position of last node along X [m]. Either this (along with ymax and zmax) or upper_bound must be specified (not both)."
    )
    ymin: float | None = Field(
        default=None,
        description="Float. Position of first node along Y [m]. Either this (along with xmin and zmin) or lower_bound must be specified (not both)."
    )
    ymax: float | None = Field(
        default=None,
        description="Float. Position of last node along Y [m]. Either this (along with xmax and zmax) or upper_bound must be specified (not both)."
    )
    zmin: float | None = Field(
        default=None,
        description="Float. Position of first node along Z [m]. Either this (along with xmin and ymin) or lower_bound must be specified (not both)."
    )
    zmax: float | None = Field(
        default=None,
        description="Float. Position of last node along Z [m]. Either this (along with xmax and ymax) or upper_bound must be specified (not both)."
    )
    bc_xmin: str | None = Field(
        default=None,
        description="String. Boundary condition at min X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_ymin and bc_zmin) or lower_boundary_conditions must be specified (not both)."
    )
    bc_xmax: str | None = Field(
        default=None,
        description="String. Boundary condition at max X: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_ymax and bc_zmax) or upper_boundary_conditions must be specified (not both)."
    )
    bc_ymin: str | None = Field(
        default=None,
        description="String. Boundary condition at min Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_xmin and bc_zmin) or lower_boundary_conditions must be specified (not both)."
    )
    bc_ymax: str | None = Field(
        default=None,
        description="String. Boundary condition at max Y: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_xmax and bc_zmax) or upper_boundary_conditions must be specified (not both)."
    )
    bc_zmin: str | None = Field(
        default=None,
        description="String. Boundary condition at min Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_xmin and bc_ymin) or lower_boundary_conditions must be specified (not both)."
    )
    bc_zmax: str | None = Field(
        default=None,
        description="String. Boundary condition at max Z: One of periodic, open, dirichlet, absorbing_silver_mueller, or neumann. Either this (along with bc_xmax and bc_ymax) or upper_boundary_conditions must be specified (not both)."
    )

    # Particle boundary parameters (vector form)
    lower_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle lower bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    upper_bound_particles: Sequence[float] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of floats, optional. Position of particle upper bound [m]. By default, if not specified, particle boundary values are the same as field boundary values."
    )
    lower_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at lower boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )
    upper_boundary_conditions_particles: Sequence[str] | None = Field(
        default=None,
        min_length=2,
        max_length=2,
        description="Vector of strings, optional. Conditions at upper boundaries for particles: periodic, absorbing, reflect or thermal. By default, if not specified, particle boundary conditions are the same as field boundary conditions."
    )

    # Particle boundary parameters (individual form)
    xmin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along X [m]."
    )
    xmax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along X [m]."
    )
    ymin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along Y [m]."
    )
    ymax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along Y [m]."
    )
    zmin_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of min particle boundary along Z [m]."
    )
    zmax_particles: float | None = Field(
        default=None,
        description="Float, optional. Position of max particle boundary along Z [m]."
    )
    bc_xmin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_xmax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_ymin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_ymax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_zmin_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal."
    )
    bc_zmax_particles: str | None = Field(
        default=None,
        description="String, optional. Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal."
    )

    # Other parameters
    moving_window_velocity: Sequence[float] | None = Field(
        default=None,
        description="Vector of floats, optional. Moving frame velocity [m/s]"
    )
    refined_regions: list[list] = Field(
        default_factory=list,
        description="List of lists, optional. List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor], with level being the refinement level (1 being first, 2 being second etc), lo and hi being vectors of length 3 specifying the extent of the region, and refinement_factor defaulting to [2,2,2] (relative to next lower level)."
    )
    guard_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of guard cells used along each direction"
    )
    pml_cells: Sequence[int] | None = Field(
        default=None,
        description="Vector of integers, optional. Number of Perfectly Matched Layer (PML) cells along each direction"
    )

    @model_validator(mode='before')
    @classmethod
    def _normalize_parameters(cls, data: Any) -> Any:
        """Normalize individual parameters to vector form before validation"""
        if not isinstance(data, dict):
            return data
        
        data = dict(data)
        
        # Validate either/or constraints for grid parameters
        number_of_cells = data.get('number_of_cells')
        nx = data.get('nx')
        ny = data.get('ny')
        nz = data.get('nz')
        if (number_of_cells is None) == ((nx is not None) and (ny is not None) and (nz is not None)):
            raise ValueError('Exactly one of number_of_cells or (nx, ny, and nz) must be specified')
        
        lower_bound = data.get('lower_bound')
        xmin = data.get('xmin')
        ymin = data.get('ymin')
        zmin = data.get('zmin')
        if (lower_bound is None) == ((xmin is not None) and (ymin is not None) and (zmin is not None)):
            raise ValueError('Exactly one of lower_bound or (xmin, ymin, and zmin) must be specified')
        
        upper_bound = data.get('upper_bound')
        xmax = data.get('xmax')
        ymax = data.get('ymax')
        zmax = data.get('zmax')
        if (upper_bound is None) == ((xmax is not None) and (ymax is not None) and (zmax is not None)):
            raise ValueError('Exactly one of upper_bound or (xmax, ymax, and zmax) must be specified')
        
        lower_boundary_conditions = data.get('lower_boundary_conditions')
        bc_xmin = data.get('bc_xmin')
        bc_ymin = data.get('bc_ymin')
        bc_zmin = data.get('bc_zmin')
        if (lower_boundary_conditions is None) == ((bc_xmin is not None) and (bc_ymin is not None) and (bc_zmin is not None)):
            raise ValueError('Exactly one of lower_boundary_conditions or (bc_xmin, bc_ymin, and bc_zmin) must be specified')
        
        upper_boundary_conditions = data.get('upper_boundary_conditions')
        bc_xmax = data.get('bc_xmax')
        bc_ymax = data.get('bc_ymax')
        bc_zmax = data.get('bc_zmax')
        if (upper_boundary_conditions is None) == ((bc_xmax is not None) and (bc_ymax is not None) and (bc_zmax is not None)):
            raise ValueError('Exactly one of upper_boundary_conditions or (bc_xmax, bc_ymax, and bc_zmax) must be specified')
        
        # Convert individual parameters to vectors if vectors not provided
        if number_of_cells is None and nx is not None and ny is not None and nz is not None:
            data['number_of_cells'] = [nx, ny, nz]
            data.pop('nx', None)
            data.pop('ny', None)
            data.pop('nz', None)
        elif number_of_cells is not None:
            data.pop('nx', None)
            data.pop('ny', None)
            data.pop('nz', None)
        
        if lower_bound is None and xmin is not None and ymin is not None and zmin is not None:
            data['lower_bound'] = [xmin, ymin, zmin]
            data.pop('xmin', None)
            data.pop('ymin', None)
            data.pop('zmin', None)
        elif lower_bound is not None:
            data.pop('xmin', None)
            data.pop('ymin', None)
            data.pop('zmin', None)
        
        if upper_bound is None and xmax is not None and ymax is not None and zmax is not None:
            data['upper_bound'] = [xmax, ymax, zmax]
            data.pop('xmax', None)
            data.pop('ymax', None)
            data.pop('zmax', None)
        elif upper_bound is not None:
            data.pop('xmax', None)
            data.pop('ymax', None)
            data.pop('zmax', None)
        
        if lower_boundary_conditions is None and bc_xmin is not None and bc_ymin is not None and bc_zmin is not None:
            data['lower_boundary_conditions'] = [bc_xmin, bc_ymin, bc_zmin]
            data.pop('bc_xmin', None)
            data.pop('bc_ymin', None)
            data.pop('bc_zmin', None)
        elif lower_boundary_conditions is not None:
            data.pop('bc_xmin', None)
            data.pop('bc_ymin', None)
            data.pop('bc_zmin', None)
        
        if upper_boundary_conditions is None and bc_xmax is not None and bc_ymax is not None and bc_zmax is not None:
            data['upper_boundary_conditions'] = [bc_xmax, bc_ymax, bc_zmax]
            data.pop('bc_xmax', None)
            data.pop('bc_ymax', None)
            data.pop('bc_zmax', None)
        elif upper_boundary_conditions is not None:
            data.pop('bc_xmax', None)
            data.pop('bc_ymax', None)
            data.pop('bc_zmax', None)
        
        # Handle particle boundaries
        if data.get('lower_bound_particles') is None:
            xmin_p = data.get('xmin_particles')
            ymin_p = data.get('ymin_particles')
            zmin_p = data.get('zmin_particles')
            if xmin_p is not None or ymin_p is not None or zmin_p is not None:
                data['lower_bound_particles'] = [xmin_p, ymin_p, zmin_p]
                data.pop('xmin_particles', None)
                data.pop('ymin_particles', None)
                data.pop('zmin_particles', None)
        if data.get('upper_bound_particles') is None:
            xmax_p = data.get('xmax_particles')
            ymax_p = data.get('ymax_particles')
            zmax_p = data.get('zmax_particles')
            if xmax_p is not None or ymax_p is not None or zmax_p is not None:
                data['upper_bound_particles'] = [xmax_p, ymax_p, zmax_p]
                data.pop('xmax_particles', None)
                data.pop('ymax_particles', None)
                data.pop('zmax_particles', None)
        if data.get('lower_boundary_conditions_particles') is None:
            bc_xmin_p = data.get('bc_xmin_particles')
            bc_ymin_p = data.get('bc_ymin_particles')
            bc_zmin_p = data.get('bc_zmin_particles')
            if bc_xmin_p is not None or bc_ymin_p is not None or bc_zmin_p is not None:
                data['lower_boundary_conditions_particles'] = [bc_xmin_p, bc_ymin_p, bc_zmin_p]
                data.pop('bc_xmin_particles', None)
                data.pop('bc_ymin_particles', None)
                data.pop('bc_zmin_particles', None)
        if data.get('upper_boundary_conditions_particles') is None:
            bc_xmax_p = data.get('bc_xmax_particles')
            bc_ymax_p = data.get('bc_ymax_particles')
            bc_zmax_p = data.get('bc_zmax_particles')
            if bc_xmax_p is not None or bc_ymax_p is not None or bc_zmax_p is not None:
                data['upper_boundary_conditions_particles'] = [bc_xmax_p, bc_ymax_p, bc_zmax_p]
                data.pop('bc_xmax_particles', None)
                data.pop('bc_ymax_particles', None)
                data.pop('bc_zmax_particles', None)
        
        return data
    
    @model_validator(mode='after')
    def _validate_and_normalize(self) -> Self:
        """Validate dimensions and normalize particle boundaries"""
        # At this point, vectors should already be set by mode='before' validator
        # Handle particle boundaries - default to field boundaries if not specified
        if self.lower_bound_particles is None:
            if self.xmin_particles is None and self.ymin_particles is None and self.zmin_particles is None:
                self.lower_bound_particles = list(self.lower_bound)
            else:
                self.lower_bound_particles = [self.xmin_particles, self.ymin_particles, self.zmin_particles]
                self.xmin_particles = None
                self.ymin_particles = None
                self.zmin_particles = None
        if self.upper_bound_particles is None:
            if self.xmax_particles is None and self.ymax_particles is None and self.zmax_particles is None:
                self.upper_bound_particles = list(self.upper_bound)
            else:
                self.upper_bound_particles = [self.xmax_particles, self.ymax_particles, self.zmax_particles]
                self.xmax_particles = None
                self.ymax_particles = None
                self.zmax_particles = None
        if self.lower_boundary_conditions_particles is None:
            if self.bc_xmin_particles is None and self.bc_ymin_particles is None and self.bc_zmin_particles is None:
                self.lower_boundary_conditions_particles = list(self.lower_boundary_conditions)
            else:
                self.lower_boundary_conditions_particles = [self.bc_xmin_particles, self.bc_ymin_particles, self.bc_zmin_particles]
                self.bc_xmin_particles = None
                self.bc_ymin_particles = None
                self.bc_zmin_particles = None
        if self.upper_boundary_conditions_particles is None:
            if self.bc_xmax_particles is None and self.bc_ymax_particles is None and self.bc_zmax_particles is None:
                self.upper_boundary_conditions_particles = list(self.upper_boundary_conditions)
            else:
                self.upper_boundary_conditions_particles = [self.bc_xmax_particles, self.bc_ymax_particles, self.bc_zmax_particles]
                self.bc_xmax_particles = None
                self.bc_ymax_particles = None
                self.bc_zmax_particles = None
        
        # Dimensions are validated by Field(min_length=3, max_length=3) constraints
        
        # Process refined regions
        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2, 2, 2])
            if len(region[1]) != 3:
                raise ValueError('The lo extent of the refined region must be a vector of length 3')
            if len(region[2]) != 3:
                raise ValueError('The hi extent of the refined region must be a vector of length 3')
            if len(region[3]) != 3:
                raise ValueError('The refinement factor of the refined region must be a vector of length 3')
        
        return self

    def add_refined_region(self, level, lo, hi, refinement_factor=[2, 2, 2]):
        """Add a refined region.

        Parameters
        ----------
        level: integer
            The refinement level, with 1 being the first level of refinement, 2 being the second etc.

        lo, hi: vectors of floats
            Each is a vector of length 3 specifying the extent of the region

        refinement_factor: vector of integers, optional
            Defaulting to [2,2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])
