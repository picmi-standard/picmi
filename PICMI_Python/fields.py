"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys
import typing
from collections.abc import Sequence

from autoclass import autoargs
from typeguard import typechecked

from .base import _ClassWithInit
from . import picmi_types

# ---------------
# Physics objects
# ---------------

@typechecked
class PICMI_ElectromagneticSolver(_ClassWithInit):
    """
    Electromagnetic field solver
      - grid: grid object to be used by the solver (grid object)
      - method: One of 'Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', or 'ECT' (string)
        - 'Yee': standard solver using the staggered Yee grid (https://doi.org/10.1109/TAP.1966.1138693)
        - 'CKC': solver with the extended Cole-Karkkainen-Cowan stencil with better dispersion properties
          (https://doi.org/10.1103/PhysRevSTAB.16.041303)
        - 'Lehe': CKC-style solver with modified dispersion (https://doi.org/10.1103/PhysRevSTAB.16.021301)
        - 'PSTD': Spectral solver with finite difference in time domain, e.g., Q. H. Liu, Letters 15 (3) (1997) 158–165
        - 'PSATD': Spectral solver with analytic in time domain (https://doi.org/10.1016/j.jcp.2013.03.010)
        - 'DS': Directional Splitting after Yasuhiko Sentoku (https://doi.org/10.1140/epjd/e2014-50162-y)
        - 'ECT': Enlarged Cell Technique solver, allowing internal conductors (https://doi.org/10.1109/APS.2005.1551259)
      - stencil_order: Order of stencil for each axis (-1=infinite) (vector of integers)
      - cfl = None: Fraction of the Courant-Friedrich-Lewy criteria [1] (float)
      - l_nodal = None: Quantities are at nodes if True, staggered otherwise (boolean)
      - source_smoother = None: Smoother object to apply to the sources (smoother object)
      - field_smoother = None: Smoother object to apply to the fields (smoother object)
      - subcycling = None: level of subcycling for the GPSTD solver (integer)
      - galilean_velocity = None: Velocity of Galilean reference frame (vector of floats) [m/s]
      - divE_cleaning = None: Solver uses div(E) cleaning if True (boolean)
      - divB_cleaning = None: Solver uses div(B) cleaning if True (boolean)
      - pml_divE_cleaning = None: Solver uses div(E) cleaning in the PML if True (boolean)
      - pml_divB_cleaning = None: Solver uses div(B) cleaning in the PML if True (boolean)
    """

    methods_list = ['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', 'ECT']

    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.FieldsGridType,
                       method : str = None,
                       stencil_order : picmi_types.VectorInt3 = None,
                       cfl : float = None,
                       l_nodal : bool = None,
                       source_smoother : picmi_types.SmootherType = None,
                       field_smoother : picmi_types.SmootherType = None,
                       subcycling : int = None,
                       galilean_velocity : picmi_types.VectorFloat3 = None,
                       divE_cleaning : bool = None,
                       divB_cleaning : bool = None,
                       pml_divE_cleaning : bool = None,
                       pml_divB_cleaning : bool = None,
                       **kw):

        assert method is None or method in PICMI_ElectromagneticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectromagneticSolver.methods_list))

        self.handle_init(kw)


@typechecked
class PICMI_ElectrostaticSolver(_ClassWithInit):
    """
    Electrostatic field solver
      - grid: grid object to be used by the solver (grid object)
      - method = None: One of 'FFT', or 'Multigrid' (string)
      - required_precision: Level of precision required for iterative solvers (float)
      - maximum_iterations: Maximum number of iterations for iterative solvers (integer)
    """

    methods_list = ['FFT', 'Multigrid']

    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.FieldsGridType,
                       method : str = None,
                       required_precision : float = None,
                       maximum_iterations : int = None,
                       **kw):

        assert method is None or method in PICMI_ElectrostaticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectrostaticSolver.methods_list))

        self.handle_init(kw)


@typechecked
class PICMI_MagnetostaticSolver(_ClassWithInit):
    """
    Magnetostatic field solver
      - grid: grid object to be used by the solver (grid object)
      - method = None: One of 'FFT', or 'Multigrid' (string)
    """

    methods_list = ['FFT', 'Multigrid']

    @autoargs(exclude=['kw'])
    def __init__(self, grid : picmi_types.FieldsGridType,
                       method : str = None,
                       **kw):

        assert method is None or method in PICMI_MagnetostaticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_MagnetostaticSolver.methods_list))

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


@typechecked
class PICMI_BinomialSmoother(_ClassWithInit):
    """
    Descibes a binomial smoother operator (applied to grids)
    - n_pass: Number of passes along each axis (vector)
    - compensation: Flags whether to apply comensation along each axis (vector of booleans)
    - stride: Stride along each axis. (vector)
    - alpha: Smoothing coefficients along each axis. (vector)
    """
    @autoargs(exclude=['kw'])
    def __init__(self, n_pass : picmi_types.VectorInt3 = None,
                       compensation : picmi_types.VectorBool3 = None,
                       stride : picmi_types.VectorInt3 = None,
                       alpha : picmi_types.VectorFloat3 = None,
                       **kw):

        self.handle_init(kw)


@typechecked
class PICMI_Cartesian1DGrid(_ClassWithInit):
    """
    One-dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

      - number_of_cells: Number of cells along each axis (number of nodes is number_of_cells+1) (vector)
      - lower_bound: Position of the node at the lower bound (vector) [m]
      - upper_bound: Position of the node at the upper bound (vector) [m]
      - lower_boundary_conditions: Conditions at lower boundaries, periodic, open, dirichlet, or neumann (vector)
      - upper_boundary_conditions: Conditions at upper boundaries, periodic, open, dirichlet, or neumann (vector)

      - nx: Number of cells along X (number of nodes=nx+1)
      - xmin: Position of first node along X [m]
      - xmax: Position of last node along X [m]
      - bc_xmin: Boundary condition at min X: One of periodic, open, dirichlet, or neumann
      - bc_xmax: Boundary condition at max X: One of periodic, open, dirichlet, or neumann

      - moving_window_velocity: Moving frame velocity (vector) [m/s]

      - refined_regions: List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
                         with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
                         lo and hi being vectors of length 1 specifying the extent of the region,
                         and refinement_factor defaulting to [2] (relative to next lower level)

      - lower_bound_particles: Position of particle lower bound (vector of floats) [m]
      - upper_bound_particles: Position of particle upper bound (vector of floats) [m]
      - xmin_particles: Position of min particle boundary along X [m] (float)
      - xmax_particles: Position of max particle boundary along X [m] (float)
      - lower_boundary_conditions_particles: Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - upper_boundary_conditions_particles: Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - bc_xmin_particles: Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_xmax_particles: Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal (string)
      - guard_cells = None: number of guard cells used along each direction (vector of integers)
      - pml_cells = None: number of Perfectly Matched Layer (PML) cells along each direction (vector of integers)
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 1

    @autoargs(exclude=['kw'])
    def __init__(self, number_of_cells : picmi_types.VectorInt1 = None,
                       lower_bound : picmi_types.VectorFloat1 = None,
                       upper_bound : picmi_types.VectorFloat1 = None,
                       lower_boundary_conditions : picmi_types.VectorString1 = None,
                       upper_boundary_conditions : picmi_types.VectorString1 = None,
                       nx : int = None,
                       xmin : float = None,
                       xmax : float = None,
                       bc_xmin : str = None,
                       bc_xmax : str = None,
                       moving_window_velocity : picmi_types.VectorFloat1 = None,
                       refined_regions : Sequence[list] = None,
                       lower_bound_particles : picmi_types.VectorFloat1 = None,
                       upper_bound_particles : picmi_types.VectorFloat1 = None,
                       xmin_particles : float = None,
                       xmax_particles : float = None,
                       lower_boundary_conditions_particles : picmi_types.VectorString1 = None,
                       upper_boundary_conditions_particles : picmi_types.VectorString1 = None,
                       bc_xmin_particles : str = None,
                       bc_xmax_particles : str = None,
                       guard_cells : picmi_types.VectorInt1 = None,
                       pml_cells : picmi_types.VectorInt1 = None,
                       **kw):

        # Sanity check and init of input arguments related to grid parameters
        assert (number_of_cells is None) and (nx is not None) or \
               (number_of_cells is not None) and (nx is None), \
                Exception('Either number_of_cells or nx must be specified')
        assert (lower_bound is None) and (xmin is not None) or \
               (lower_bound is not None) and (xmin is None), \
                Exception('Either lower_bound or xmin must be specified')
        assert (upper_bound is None) and (xmax is not None) or \
               (upper_bound is not None) and (xmax is None), \
                Exception('Either upper_bound or xmax must be specified')
        assert (lower_boundary_conditions is None) and (bc_xmin is not None) or \
               (lower_boundary_conditions is not None) and (bc_xmin is None), \
                Exception('Either lower_boundary_conditions or bc_xmin')
        assert (upper_boundary_conditions is None) and (bc_xmax is not None) or \
               (upper_boundary_conditions is not None) and (bc_xmax is None), \
                Exception('Either upper_boundary_conditions or bc_xmax must be specified')

        if number_of_cells is None:
            self.number_of_cells = [nx]
        if lower_bound is None:
            self.lower_bound = [xmin]
        if upper_bound is None:
            self.upper_bound = [xmax]
        if lower_boundary_conditions is None:
            self.lower_boundary_conditions = [bc_xmin]
        if upper_boundary_conditions is None:
            self.upper_boundary_conditions = [bc_xmax,]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (xmin_particles is None):
                self.lower_bound_particles = lower_bound
            else:
                self.lower_bound_particles = [xmin_particles]
        if upper_bound_particles is None:
            if (xmax_particles is None):
                self.upper_bound_particles = upper_bound
            else:
                self.upper_bound_particles = [xmax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_xmin_particles is None):
                self.lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [bc_xmin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_xmax_particles is None):
                self.upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [bc_xmax_particles]

        self.refined_regions = []
        if refined_regions is not None:
            for region in refined_regions:
                self.add_refined_region(*region)

        self.handle_init(kw)

    def add_refined_region(self, level : int,
                                 lo : picmi_types.VectorInt1,
                                 hi : picmi_types.VectorInt1,
                                 refinement_factor : picmi_types.VectorInt1 = None):
        """Add a refined region.
        - level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        - lo, hi: vectors of length 1 specifying the extent of the region
        - refinement_factor: defaulting to [2] (relative to next lower level)
        """
        if refinement_factor is None:
            refinement_factor = [2]
        assert len(lo) == 1, Exception('The lo extent of the refined region must be a vector of length 1')
        assert len(hi) == 1, Exception('The hi extent of the refined region must be a vector of length 1')
        assert len(refinement_factor) == 1, Exception('The refinement factor of the refined region must be a vector of length 1')
        self.refined_regions.append([level, lo, hi, refinement_factor])


@typechecked
class PICMI_CylindricalGrid(_ClassWithInit):
    """
    Axisymmetric, cylindrical grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

      - number_of_cells: Number of cells along each axis (number of nodes is number_of_cells+1) (vector)
      - lower_bound: Position of the node at the lower bound (vector) [m]
      - upper_bound: Position of the node at the upper bound (vector) [m]
      - lower_boundary_conditions: Conditions at lower boundaries, periodic, open, dirichlet, or neumann (vector)
      - upper_boundary_conditions: Conditions at upper boundaries, periodic, open, dirichlet, or neumann (vector)

      - nr: Number of cells along R (number of nodes=nr+1)
      - nz: Number of cells along Z (number of nodes=nz+1)
      - n_azimuthal_modes: Number of azimuthal modes
      - rmin: Position of first node along R [m]
      - rmax: Position of last node along R [m]
      - zmin: Position of first node along Z [m]
      - zmax: Position of last node along Z [m]
      - bc_rmin: Boundary condition at min R: One of open, dirichlet, or neumann
      - bc_rmax: Boundary condition at max R: One of open, dirichlet, or neumann
      - bc_zmin: Boundary condition at min Z: One of periodic, open, dirichlet, or neumann
      - bc_zmax: Boundary condition at max Z: One of periodic, open, dirichlet, or neumann

      - moving_window_velocity: Moving frame velocity (vector) [m/s]

      - refined_regions: List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
                         with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
                         lo and hi being vectors of length 2 specifying the extent of the region,
                         and refinement_factor defaulting to [2,2] (relative to next lower level)

      - lower_bound_particles: Position of particle lower bound (vector of floats) [m]
      - upper_bound_particles: Position of particle upper bound (vector of floats) [m]
      - rmin_particles: Position of min particle boundary along R [m] (float)
      - rmax_particles: Position of max particle boundary along R [m] (float)
      - zmin_particles: Position of min particle boundary along Z [m] (float)
      - zmax_particles: Position of max particle boundary along Z [m] (float)
      - lower_boundary_conditions_particles: Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - upper_boundary_conditions_particles: Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - bc_rmin_particles: Boundary condition at min R for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_rmax_particles: Boundary condition at max R for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_zmin_particles: Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_zmax_particles: Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal (string)
      - guard_cells = None: number of guard cells used along each direction (vector of integers)
      - pml_cells = None: number of Perfectly Matched Layer (PML) cells along each direction (vector of integers)
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nr and nz) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 2

    @autoargs(exclude=['kw'])
    def __init__(self, number_of_cells : picmi_types.VectorInt2 = None,
                       lower_bound : picmi_types.VectorFloat2 = None,
                       upper_bound : picmi_types.VectorFloat2 = None,
                       lower_boundary_conditions : picmi_types.VectorString2 = None,
                       upper_boundary_conditions : picmi_types.VectorString2 = None,
                       nr : int = None,
                       nz : int = None,
                       n_azimuthal_modes : int = None,
                       rmin : float = None,
                       rmax : float = None,
                       zmin : float = None,
                       zmax : float = None,
                       bc_rmin : str = None,
                       bc_rmax : str = None,
                       bc_zmin : str = None,
                       bc_zmax : str = None,
                       moving_window_velocity : picmi_types.VectorFloat2 = None,
                       refined_regions : Sequence[list] = None,
                       lower_bound_particles : picmi_types.VectorFloat2 = None,
                       upper_bound_particles : picmi_types.VectorFloat2 = None,
                       rmin_particles : float = None,
                       rmax_particles : float = None,
                       zmin_particles : float = None,
                       zmax_particles : float = None,
                       lower_boundary_conditions_particles : picmi_types.VectorString2 = None,
                       upper_boundary_conditions_particles : picmi_types.VectorString2 = None,
                       bc_rmin_particles : str = None,
                       bc_rmax_particles : str = None,
                       bc_zmin_particles : str = None,
                       bc_zmax_particles : str = None,
                       guard_cells : picmi_types.VectorInt2 = None,
                       pml_cells : picmi_types.VectorInt2 = None,
                       **kw):

        # Sanity check and init of input arguments related to grid parameters
        assert (number_of_cells is None) and (nr is not None and nz is not None) or \
               (number_of_cells is not None) and (nr is None and nz is None), \
                Exception('Either number_of_cells or nr and nz must be specified')
        assert (lower_bound is None) and (rmin is not None and zmin is not None) or \
               (lower_bound is not None) and (rmin is None and zmin is None), \
                Exception('Either lower_bound or rmin and zmin must be specified')
        assert (upper_bound is None) and (rmax is not None and zmax is not None) or \
               (upper_bound is not None) and (rmax is None and zmax is None), \
                Exception('Either upper_bound or rmax and zmax must be specified')
        # --- Allow bc_rmin to be None since it will usually be the axis.
        assert (lower_boundary_conditions is None) and (bc_zmin is not None) or \
               (lower_boundary_conditions is not None) and (bc_rmin is None and bc_zmin is None), \
                Exception('Either lower_boundary_conditions or bc_rmin and bc_zmin must be specified')
        assert (upper_boundary_conditions is None) and (bc_rmax is not None and bc_zmax is not None) or \
               (upper_boundary_conditions is not None) and (bc_rmax is None and bc_zmax is None), \
                Exception('Either upper_boundary_conditions or bc_rmax and bc_zmax must be specified')

        if number_of_cells is None:
            self.number_of_cells = [nr, nz]
        if lower_bound is None:
            self.lower_bound = [rmin, zmin]
        if upper_bound is None:
            self.upper_bound = [rmax, zmax]
        if lower_boundary_conditions is None:
            self.lower_boundary_conditions = [bc_rmin, bc_zmin]
        if upper_boundary_conditions is None:
            self.upper_boundary_conditions = [bc_rmax, bc_zmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (rmin_particles is None) and (zmin_particles is None):
                self.lower_bound_particles = lower_bound
            else:
                self.lower_bound_particles = [rmin_particles, zmin_particles]
        if upper_bound_particles is None:
            if (rmax_particles is None) and (zmax_particles is None):
                self.upper_bound_particles = upper_bound
            else:
                self.upper_bound_particles = [rmax_particles, zmax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_rmin_particles is None) and (bc_zmin_particles is None):
                self.lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [bc_rmin_particles, bc_zmin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_rmax_particles is None) and (bc_zmax_particles is None):
                self.upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [bc_rmax_particles, bc_zmax_particles]

        self.refined_regions = []
        if refined_regions is not None:
            for region in refined_regions:
                self.add_refined_region(*region)

        self.handle_init(kw)

    def add_refined_region(self, level : int,
                                 lo : picmi_types.VectorInt2,
                                 hi : picmi_types.VectorInt2,
                                 refinement_factor : picmi_types.VectorInt2 = None):
        """Add a refined region.
        - level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        - lo, hi: vectors of length 2 specifying the extent of the region
        - refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        if refinement_factor is None:
            refinement_factor = [2,2]
        assert len(lo) == 2, Exception('The lo extent of the refined region must be a vector of length 2')
        assert len(hi) == 2, Exception('The hi extent of the refined region must be a vector of length 2')
        assert len(refinement_factor) == 2, Exception('The refinement factor of the refined region must be a vector of length 2')
        self.refined_regions.append([level, lo, hi, refinement_factor])


@typechecked
class PICMI_Cartesian2DGrid(_ClassWithInit):
    """
    Two dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

      - number_of_cells: Number of cells along each axis (number of nodes is number_of_cells+1) (vector)
      - lower_bound: Position of the node at the lower bound (vector) [m]
      - upper_bound: Position of the node at the upper bound (vector) [m]
      - lower_boundary_conditions: Conditions at lower boundaries, periodic, open, dirichlet, or neumann (vector)
      - upper_boundary_conditions: Conditions at upper boundaries, periodic, open, dirichlet, or neumann (vector)

      - nx: Number of cells along X (number of nodes=nx+1)
      - ny: Number of cells along Y (number of nodes=ny+1)
      - xmin: Position of first node along X [m]
      - xmax: Position of last node along X [m]
      - ymin: Position of first node along Y [m]
      - ymax: Position of last node along Y [m]
      - bc_xmin: Boundary condition at min X: One of periodic, open, dirichlet, or neumann
      - bc_xmax: Boundary condition at max X: One of periodic, open, dirichlet, or neumann
      - bc_ymin: Boundary condition at min Y: One of periodic, open, dirichlet, or neumann
      - bc_ymax: Boundary condition at max Y: One of periodic, open, dirichlet, or neumann

      - moving_window_velocity: Moving frame velocity (vector) [m/s]

      - refined_regions: List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
                         with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
                         lo and hi being vectors of length 2 specifying the extent of the region,
                         and refinement_factor defaulting to [2,2] (relative to next lower level)

      - lower_bound_particles: Position of particle lower bound (vector of floats) [m]
      - upper_bound_particles: Position of particle upper bound (vector of floats) [m]
      - xmin_particles: Position of min particle boundary along X [m] (float)
      - xmax_particles: Position of max particle boundary along X [m] (float)
      - ymin_particles: Position of min particle boundary along Y [m] (float)
      - ymax_particles: Position of max particle boundary along Y [m] (float)
      - lower_boundary_conditions_particles: Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - upper_boundary_conditions_particles: Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - bc_xmin_particles: Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_xmax_particles: Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_ymin_particles: Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_ymax_particles: Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal (string)
      - guard_cells = None: number of guard cells used along each direction (vector of integers)
      - pml_cells = None: number of Perfectly Matched Layer (PML) cells along each direction (vector of integers)
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx and ny) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 2

    @autoargs(exclude=['kw'])
    def __init__(self, number_of_cells : picmi_types.VectorInt2 = None,
                       lower_bound : picmi_types.VectorFloat2 = None,
                       upper_bound : picmi_types.VectorFloat2 = None,
                       lower_boundary_conditions : picmi_types.VectorString2 = None,
                       upper_boundary_conditions : picmi_types.VectorString2 = None,
                       nx : int = None,
                       ny : int = None,
                       xmin : float = None,
                       xmax : float = None,
                       ymin : float = None,
                       ymax : float = None,
                       bc_xmin : str = None,
                       bc_xmax : str = None,
                       bc_ymin : str = None,
                       bc_ymax : str = None,
                       moving_window_velocity : picmi_types.VectorFloat2 = None,
                       refined_regions : Sequence[list] = None,
                       lower_bound_particles : picmi_types.VectorFloat2 = None,
                       upper_bound_particles : picmi_types.VectorFloat2 = None,
                       xmin_particles : float = None,
                       xmax_particles : float = None,
                       ymin_particles : float = None,
                       ymax_particles : float = None,
                       lower_boundary_conditions_particles : picmi_types.VectorString2 = None,
                       upper_boundary_conditions_particles : picmi_types.VectorString2 = None,
                       bc_xmin_particles : str = None,
                       bc_xmax_particles : str = None,
                       bc_ymin_particles : str = None,
                       bc_ymax_particles : str = None,
                       guard_cells : picmi_types.VectorInt2 = None,
                       pml_cells : picmi_types.VectorInt2 = None,
                       **kw):

        # Sanity check and init of input arguments related to grid parameters
        assert (number_of_cells is None) and (nx is not None and ny is not None) or \
               (number_of_cells is not None) and (nx is None and ny is None), \
                Exception('Either number_of_cells or nx and ny must be specified')
        assert (lower_bound is None) and (xmin is not None and ymin is not None) or \
               (lower_bound is not None) and (xmin is None and ymin is None), \
                Exception('Either lower_bound or xmin and ymin must be specified')
        assert (upper_bound is None) and (xmax is not None and ymax is not None) or \
               (upper_bound is not None) and (xmax is None and ymax is None), \
                Exception('Either upper_bound or xmax and ymax must be specified')
        assert (lower_boundary_conditions is None) and (bc_xmin is not None and bc_ymin is not None) or \
               (lower_boundary_conditions is not None) and (bc_xmin is None and bc_ymin is None), \
                Exception('Either lower_boundary_conditions or bc_xmin and bc_ymin must be specified')
        assert (upper_boundary_conditions is None) and (bc_xmax is not None and bc_ymax is not None) or \
               (upper_boundary_conditions is not None) and (bc_xmax is None and bc_ymax is None), \
                Exception('Either upper_boundary_conditions or bc_xmax and bc_ymax must be specified')

        if number_of_cells is None:
            self.number_of_cells = [nx, ny]
        if lower_bound is None:
            self.lower_bound = [xmin, ymin]
        if upper_bound is None:
            self.upper_bound = [xmax, ymax]
        if lower_boundary_conditions is None:
            self.lower_boundary_conditions = [bc_xmin, bc_ymin]
        if upper_boundary_conditions is None:
            self.upper_boundary_conditions = [bc_xmax, bc_ymax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (xmin_particles is None) and (ymin_particles is None):
                self.lower_bound_particles = lower_bound
            else:
                self.lower_bound_particles = [xmin_particles, ymin_particles]
        if upper_bound_particles is None:
            if (xmax_particles is None) and (ymax_particles is None):
                self.upper_bound_particles = upper_bound
            else:
                self.upper_bound_particles = [xmax_particles, ymax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_xmin_particles is None) and (bc_ymin_particles is None):
                self.lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [bc_xmin_particles, bc_ymin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_xmax_particles is None) and (bc_ymax_particles is None):
                self.upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [bc_xmax_particles, bc_ymax_particles]

        self.refined_regions = []
        if refined_regions is not None:
            for region in refined_regions:
                self.add_refined_region(*region)

        self.handle_init(kw)

    def add_refined_region(self, level : int,
                                 lo : picmi_types.VectorInt2,
                                 hi : picmi_types.VectorInt2,
                                 refinement_factor : picmi_types.VectorInt2 = None):
        """Add a refined region.
        - level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        - lo, hi: vectors of length 2 specifying the extent of the region
        - refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        if refinement_factor is None:
            refinement_factor = [2,2]
        assert len(lo) == 2, Exception('The lo extent of the refined region must be a vector of length 2')
        assert len(hi) == 2, Exception('The hi extent of the refined region must be a vector of length 2')
        assert len(refinement_factor) == 2, Exception('The refinement factor of the refined region must be a vector of length 2')
        self.refined_regions.append([level, lo, hi, refinement_factor])


@typechecked
class PICMI_Cartesian3DGrid(_ClassWithInit):
    """
    Three dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

      - number_of_cells: Number of cells along each axis (number of nodes is number_of_cells+1) (vector of integers)
      - lower_bound: Position of the node at the lower bound (vector of floats) [m]
      - upper_bound: Position of the node at the upper bound (vector of floats) [m]
      - lower_boundary_conditions: Conditions at lower boundaries, periodic, open, dirichlet, or neumann (vector of strings)
      - upper_boundary_conditions: Conditions at upper boundaries, periodic, open, dirichlet, or neumann (vector of strings)

      - nx: Number of cells along X (number of nodes=nx+1) (integer)
      - ny: Number of cells along Y (number of nodes=ny+1) (integer)
      - nz: Number of cells along Z (number of nodes=nz+1) (integer)
      - xmin: Position of first node along X [m] (float)
      - xmax: Position of last node along X [m] (float)
      - ymin: Position of first node along Y [m] (float)
      - ymax: Position of last node along Y [m] (float)
      - zmin: Position of first node along Z [m] (float)
      - zmax: Position of last node along Z [m] (float)
      - bc_xmin: Boundary condition at min X: One of periodic, open, dirichlet, or neumann (string)
      - bc_xmax: Boundary condition at max X: One of periodic, open, dirichlet, or neumann (string)
      - bc_ymin: Boundary condition at min Y: One of periodic, open, dirichlet, or neumann (string)
      - bc_ymax: Boundary condition at max Y: One of periodic, open, dirichlet, or neumann (string)
      - bc_zmin: Boundary condition at min Z: One of periodic, open, dirichlet, or neumann (string)
      - bc_zmax: Boundary condition at max Z: One of periodic, open, dirichlet, or neumann (string)
      - moving_window_velocity: Moving frame velocity (vector) [m/s]

      - refined_regions: List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
                         with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
                         lo and hi being vectors of length 3 specifying the extent of the region,
                         and refinement_factor defaulting to [2,2,2] (relative to next lower level)

      - lower_bound_particles: Position of particle lower bound (vector of floats) [m]
      - upper_bound_particles: Position of particle upper bound (vector of floats) [m]
      - xmin_particles: Position of min particle boundary along X [m] (float)
      - xmax_particles: Position of max particle boundary along X [m] (float)
      - ymin_particles: Position of min particle boundary along Y [m] (float)
      - ymax_particles: Position of max particle boundary along Y [m] (float)
      - zmin_particles  Position of min particle boundary along Z [m] (float)
      - zmax_particles: Position of max particle boundary along Z [m] (float)
      - lower_boundary_conditions_particles: Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - upper_boundary_conditions_particles: Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal (vector of strings)
      - bc_xmin_particles: Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_xmax_particles: Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_ymin_particles: Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_ymax_particles: Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_zmin_particles: Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal (string)
      - bc_zmax_particles: Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal (string)
      - guard_cells = None: number of guard cells used along each direction (vector of integers)
      - pml_cells = None: number of Perfectly Matched Layer (PML) cells along each direction (vector of integers)
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx, ny, and nz) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 3

    @autoargs(exclude=['kw'])
    def __init__(self, number_of_cells : picmi_types.VectorInt3 = None,
                       lower_bound : picmi_types.VectorFloat3 = None,
                       upper_bound : picmi_types.VectorFloat3 = None,
                       lower_boundary_conditions : picmi_types.VectorString3 = None,
                       upper_boundary_conditions : picmi_types.VectorString3 = None,
                       nx : int = None,
                       ny : int = None,
                       nz : int = None,
                       xmin : float = None,
                       xmax : float = None,
                       ymin : float = None,
                       ymax : float = None,
                       zmin : float = None,
                       zmax : float = None,
                       bc_xmin : str = None,
                       bc_xmax : str = None,
                       bc_ymin : str = None,
                       bc_ymax : str = None,
                       bc_zmin : str = None,
                       bc_zmax : str = None,
                       moving_window_velocity : picmi_types.VectorFloat3 = None,
                       refined_regions : Sequence[list] = None,
                       lower_bound_particles : picmi_types.VectorFloat3 = None,
                       upper_bound_particles : picmi_types.VectorFloat3 = None,
                       xmin_particles : float = None,
                       xmax_particles : float = None,
                       ymin_particles : float = None,
                       ymax_particles : float = None,
                       zmin_particles : float = None,
                       zmax_particles : float = None,
                       lower_boundary_conditions_particles : picmi_types.VectorString3 = None,
                       upper_boundary_conditions_particles : picmi_types.VectorString3 = None,
                       bc_xmin_particles : str = None,
                       bc_xmax_particles : str = None,
                       bc_ymin_particles : str = None,
                       bc_ymax_particles : str = None,
                       bc_zmin_particles : str = None,
                       bc_zmax_particles : str = None,
                       guard_cells : picmi_types.VectorInt3 = None,
                       pml_cells : picmi_types.VectorInt3 = None,
                       **kw):

        # Sanity check and init of input arguments related to grid parameters
        assert (number_of_cells is None) and (nx is not None and ny is not None and nz is not None) or \
               (number_of_cells is not None) and (nx is None and ny is None and nz is None), \
                Exception('Either number_of_cells or nx, ny, and nz must be specified')
        assert (lower_bound is None) and (xmin is not None and ymin is not None and zmin is not None) or \
               (lower_bound is not None) and (xmin is None and ymin is None and zmin is None), \
                Exception('Either lower_bound or xmin, ymin, and zmin must be specified')
        assert (upper_bound is None) and (xmax is not None and ymax is not None and zmax is not None) or \
               (upper_bound is not None) and (xmax is None and ymax is None and zmax is None), \
                Exception('Either upper_bound or xmax, ymax, and zmax must be specified')
        assert (lower_boundary_conditions is None) and (bc_xmin is not None and bc_ymin is not None and bc_zmin is not None) or \
               (lower_boundary_conditions is not None) and (bc_xmin is None and bc_ymin is None and bc_zmin is None), \
                Exception('Either lower_boundary_conditions or bc_xmin, bc_ymin, and bc_zmin must be specified')
        assert (upper_boundary_conditions is None) and (bc_xmax is not None and bc_ymax is not None and bc_zmax is not None) or \
               (upper_boundary_conditions is not None) and (bc_xmax is None and bc_ymax is None and bc_zmax is None), \
                Exception('Either upper_boundary_conditions or bc_xmax, bc_ymax, and bc_zmax must be specified')

        if number_of_cells is None:
            self.number_of_cells = [nx, ny, nz]
        if lower_bound is None:
            self.lower_bound = [xmin, ymin, zmin]
        if upper_bound is None:
            self.upper_bound = [xmax, ymax, zmax]
        if lower_boundary_conditions is None:
            self.lower_boundary_conditions = [bc_xmin, bc_ymin, bc_zmin]
        if upper_boundary_conditions is None:
            self.upper_boundary_conditions = [bc_xmax, bc_ymax, bc_zmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (xmin_particles is None) and (ymin_particles is None) and (zmin_particles is None):
                self.lower_bound_particles = lower_bound
            else:
                self.lower_bound_particles = [xmin_particles, ymin_particles, zmin_particles]
        if upper_bound_particles is None:
            if (xmax_particles is None) and (ymax_particles is None) and (zmax_particles is None):
                self.upper_bound_particles = upper_bound
            else:
                self.upper_bound_particles = [xmax_particles, ymax_particles, zmax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_xmin_particles is None) and (bc_ymin_particles is None) and (bc_zmin_particles is None):
                self.lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                self.lower_boundary_conditions_particles = [bc_xmin_particles, bc_ymin_particles, bc_zmin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_xmax_particles is None) and (bc_ymax_particles is None) and (bc_zmax_particles is None):
                self.upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                self.upper_boundary_conditions_particles = [bc_xmax_particles, bc_ymax_particles, bc_zmax_particles]

        self.refined_regions = []
        if refined_regions is not None:
            for region in refined_regions:
                self.add_refined_region(*region)

        self.handle_init(kw)

    def add_refined_region(self, level : int,
                                 lo : picmi_types.VectorInt3,
                                 hi : picmi_types.VectorInt3,
                                 refinement_factor : picmi_types.VectorInt3 = None):
        """Add a refined region.
        - level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        - lo, hi: vectors of length 3 specifying the extent of the region
        - refinement_factor: defaulting to [2,2,2] (relative to next lower level)
        """
        if refinement_factor is None:
            refinement_factor = [2,2,2]
        assert len(lo) == 3, Exception('The lo extent of the refined region must be a vector of length 3')
        assert len(hi) == 3, Exception('The hi extent of the refined region must be a vector of length 3')
        assert len(refinement_factor) == 3, Exception('The refinement factor of the refined region must be a vector of length 3')
        self.refined_regions.append([level, lo, hi, refinement_factor])
