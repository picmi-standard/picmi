"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys

from .base import _ClassWithInit

# ---------------
# Physics objects
# ---------------

class PICMI_ElectromagneticSolver(_ClassWithInit):
    """
    Electromagnetic field solver

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    method: {'Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', 'ECT'}
        The advance method use to solve Maxwell's equations. The default method is code dependent.

        - 'Yee': standard solver using the staggered Yee grid (https://doi.org/10.1109/TAP.1966.1138693)

        - 'CKC': solver with the extended Cole-Karkkainen-Cowan stencil with better dispersion properties
          (https://doi.org/10.1103/PhysRevSTAB.16.041303)

        - 'Lehe': CKC-style solver with modified dispersion (https://doi.org/10.1103/PhysRevSTAB.16.021301)

        - 'PSTD': Spectral solver with finite difference in time domain, e.g., Q. H. Liu, Letters 15 (3) (1997) 158–165

        - 'PSATD': Spectral solver with analytic in time domain (https://doi.org/10.1016/j.jcp.2013.03.010)

        - 'DS': Directional Splitting after Yasuhiko Sentoku (https://doi.org/10.1140/epjd/e2014-50162-y)

        - 'ECT': Enlarged Cell Technique solver, allowing internal conductors (https://doi.org/10.1109/APS.2005.1551259)

    stencil_order: vector of integers
        Order of stencil for each axis (-1=infinite)

    cfl: float, optional
        Fraction of the Courant-Friedrich-Lewy criteria [1]

    l_nodal: bool, optional
        Quantities are at nodes if True, staggered otherwise

    source_smoother: smoother instance, optional
        Smoother object to apply to the sources

    field_smoother: smoother instance, optional
        Smoother object to apply to the fields

    subcycling: integer, optional
        Level of subcycling for the GPSTD solver

    galilean_velocity: vector of floats, optional
        Velocity of Galilean reference frame [m/s]

    divE_cleaning: bool, optional
        Solver uses div(E) cleaning if True

    divB_cleaning: bool, optional
        Solver uses div(B) cleaning if True

    pml_divE_cleaning: bool, optional
        Solver uses div(E) cleaning in the PML if True

    pml_divB_cleaning: bool, optional
        Solver uses div(B) cleaning in the PML if True
    """

    methods_list = ['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD', 'DS', 'ECT']

    def __init__(self, grid, method=None, stencil_order=None, cfl=None, l_nodal=None,
                 source_smoother=None, field_smoother=None, subcycling=None,
                 galilean_velocity=None, divE_cleaning=None, divB_cleaning=None,
                 pml_divE_cleaning=None, pml_divB_cleaning=None, **kw):

        assert method is None or method in PICMI_ElectromagneticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectromagneticSolver.methods_list))

        self.grid = grid
        self.method = method
        self.cfl = cfl
        self.stencil_order = stencil_order
        self.l_nodal = l_nodal
        self.source_smoother = source_smoother
        self.field_smoother = field_smoother
        self.subcycling = subcycling
        self.galilean_velocity = galilean_velocity
        self.divE_cleaning = divE_cleaning
        self.divB_cleaning = divB_cleaning
        self.pml_divE_cleaning = pml_divE_cleaning
        self.pml_divB_cleaning = pml_divB_cleaning

        self.handle_init(kw)


class PICMI_ElectrostaticSolver(_ClassWithInit):
    """
    Electrostatic field solver

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    method: string
        One of 'FFT', or 'Multigrid'

    required_precision: float, optional
        Level of precision required for iterative solvers

    maximum_iterations: integer, optional
        Maximum number of iterations for iterative solvers
    """

    methods_list = ['FFT', 'Multigrid']

    def __init__(self, grid, method=None,
                 required_precision=None, maximum_iterations=None, **kw):

        assert method is None or method in PICMI_ElectrostaticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectrostaticSolver.methods_list))

        self.grid = grid
        self.method = method
        self.required_precision = required_precision
        self.maximum_iterations = maximum_iterations

        self.handle_init(kw)


class PICMI_MagnetostaticSolver(_ClassWithInit):
    """
    Magnetostatic field solver

    Parameters
    ----------
    grid: grid instance
        Grid object for the diagnostic

    method: string
        One of 'FFT', or 'Multigrid'
    """

    methods_list = ['FFT', 'Multigrid']

    def __init__(self, grid, method=None, **kw):

        assert method is None or method in PICMI_MagnetostaticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_MagnetostaticSolver.methods_list))

        self.grid = grid
        self.method = method

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


class PICMI_BinomialSmoother(_ClassWithInit):
    """
    Descibes a binomial smoother operator (applied to grids)

    Parameters
    ----------
    n_pass: vector of integers
        Number of passes along each axis

    compensation: vector of booleans, optional
        Flags whether to apply comensation along each axis

    stride: vector of integers, optional
        Stride along each axis

    alpha: vector of floats, optional
        Smoothing coefficients along each axis
    """
    def __init__(self, n_pass=None, compensation=None, stride=None, alpha=None, **kw):
        self.n_pass = n_pass
        self.compensation = compensation
        self.stride = stride
        self.alpha = alpha

        self.handle_init(kw)


class PICMI_Cartesian1DGrid(_ClassWithInit):
    """
    One-dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, or neumann

    nx: integer
        Number of cells along X (number of nodes=nx+1)

    xmin: float
        Position of first node along X [m]

    xmax: float
        Position of last node along X [m]

    bc_xmin: vector of strings
        Boundary condition at min X: One of periodic, open, dirichlet, or neumann

    bc_xmax: vector of strings
        Boundary condition at max X: One of periodic, open, dirichlet, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 2 specifying the extent of the region,
        and refinement_factor defaulting to [2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    xmin_particles: float, optional
        Position of min particle boundary along X [m]

    xmax_particles: float, optional
        Position of max particle boundary along X [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_xmin_particles: string, optional
        Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal

    bc_xmax_particles: string, optional
        Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 1

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nx=None, xmin=None, xmax=None, bc_xmin=None, bc_xmax=None,
                 moving_window_velocity=None, refined_regions=[],lower_bound_particles=None, upper_bound_particles=None,
                 xmin_particles=None, xmax_particles=None,
                 lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None,
                 bc_xmin_particles=None, bc_xmax_particles=None,
                 guard_cells=None, pml_cells=None,
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
            number_of_cells = [nx]
        if lower_bound is None:
            lower_bound = [xmin]
        if upper_bound is None:
            upper_bound = [xmax]
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_xmin]
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_xmax,]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (xmin_particles is None):
                lower_bound_particles = lower_bound
            else:
                lower_bound_particles = [xmin_particles]
        if upper_bound_particles is None:
            if (xmax_particles is None):
                upper_bound_particles = upper_bound
            else:
                upper_bound_particles=[xmax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_xmin_particles is None):
                lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                lower_boundary_conditions_particles = [bc_xmin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_xmax_particles is None):
                upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                upper_boundary_conditions_particles = [bc_xmax_particles]

        # Sanity check on dimensionality of vector quantities
        assert len(number_of_cells) == 1, Exception('Wrong number of cells specified')
        assert len(lower_bound) == 1, Exception('Wrong number of lower bounds specified')
        assert len(upper_bound) == 1, Exception('Wrong number of upper bounds specified')
        assert len(lower_boundary_conditions) == 1, Exception('Wrong number of lower boundary conditions specified')
        assert len(upper_boundary_conditions) == 1, Exception('Wrong number of upper boundary conditions specified')
        assert len(lower_bound_particles) == 1, Exception('Wrong number of particle lower bounds specified')
        assert len(upper_bound_particles) == 1, Exception('Wrong number of particle upper bounds specified')
        assert len(lower_boundary_conditions_particles) == 1, Exception('Wrong number of lower particle boundary conditions specified')
        assert len(upper_boundary_conditions_particles) == 1, Exception('Wrong number of upper particle boundary conditions specified')

        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions
        self.lower_bound_particles = lower_bound_particles
        self.upper_bound_particles = upper_bound_particles
        self.lower_boundary_conditions_particles = lower_boundary_conditions_particles
        self.upper_boundary_conditions_particles = upper_boundary_conditions_particles
        self.guard_cells = guard_cells
        self.pml_cells = pml_cells

        self.moving_window_velocity = moving_window_velocity

        self.refined_regions = refined_regions

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2])
            assert len(region[1]) == 1, Exception('The lo extent of the refined region must be a vector of length 2')
            assert len(region[2]) == 1, Exception('The hi extent of the refined region must be a vector of length 2')
            assert len(region[3]) == 1, Exception('The refinement factor of the refined region must be a vector of length 2')

        self.handle_init(kw)

    def add_refined_region(self, level, lo, hi, refinement_factor=[2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_CylindricalGrid(_ClassWithInit):
    """
    Axisymmetric, cylindrical grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, or neumann

    nr: integer
        Number of cells along R (number of nodes=nr+1)

    nz: integer
        Number of cells along Z (number of nodes=nz+1)

    n_azimuthal_modes: integer
        Number of azimuthal modes

    rmin: float
        Position of first node along R [m]

    rmax: float
        Position of last node along R [m]

    zmin: float
        Position of first node along Z [m]

    zmax: float
        Position of last node along Z [m]

    bc_rmin: vector of strings
        Boundary condition at min R: One of open, dirichlet, or neumann

    bc_rmax: vector of strings
        Boundary condition at max R: One of open, dirichlet, or neumann

    bc_zmin: vector of strings
        Boundary condition at min Z: One of periodic, open, dirichlet, or neumann

    bc_zmax: vector of strings
        Boundary condition at max Z: One of periodic, open, dirichlet, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 2 specifying the extent of the region,
        and refinement_factor defaulting to [2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    rmin_particles: float, optional
        Position of min particle boundary along R [m]

    rmax_particles: float, optional
        Position of max particle boundary along R [m]

    zmin_particles: float, optional
        Position of min particle boundary along Z [m]

    zmax_particles: float, optional
        Position of max particle boundary along Z [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_rmin_particles: string, optional
        Boundary condition at min R for particles: One of periodic, absorbing, reflect, thermal

    bc_rmax_particles: string, optional
        Boundary condition at max R for particles: One of periodic, absorbing, reflect, thermal

    bc_zmin_particles: string, optional
        Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal

    bc_zmax_particles: string, optional
        Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nr and nz) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 2

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nr=None, nz=None, n_azimuthal_modes=None,
                 rmin=None, rmax=None, zmin=None, zmax=None,
                 bc_rmin=None, bc_rmax=None, bc_zmin=None, bc_zmax=None,
                 moving_window_velocity=None, refined_regions=[],
                 lower_bound_particles=None, upper_bound_particles=None,
                 rmin_particles=None, rmax_particles=None, zmin_particles=None, zmax_particles=None,
                 lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None,
                 bc_rmin_particles=None, bc_rmax_particles=None, bc_zmin_particles=None, bc_zmax_particles=None,
                 guard_cells=None, pml_cells=None, **kw):

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
        # --Allow bc_rmin to be None since it will usually be the axis.
        assert (lower_boundary_conditions is None) and (bc_zmin is not None) or \
               (lower_boundary_conditions is not None) and (bc_rmin is None and bc_zmin is None), \
                Exception('Either lower_boundary_conditions or bc_rmin and bc_zmin must be specified')
        assert (upper_boundary_conditions is None) and (bc_rmax is not None and bc_zmax is not None) or \
               (upper_boundary_conditions is not None) and (bc_rmax is None and bc_zmax is None), \
                Exception('Either upper_boundary_conditions or bc_rmax and bc_zmax must be specified')

        if number_of_cells is None:
            number_of_cells = [nr, nz]
        if lower_bound is None:
            lower_bound = [rmin, zmin]
        if upper_bound is None:
            upper_bound = [rmax, zmax]
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_rmin, bc_zmin]
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_rmax, bc_zmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (rmin_particles is None) and (zmin_particles is None):
                lower_bound_particles = lower_bound
            else:
                lower_bound_particles = [rmin_particles, zmin_particles]
        if upper_bound_particles is None:
            if (rmax_particles is None) and (zmax_particles is None):
                upper_bound_particles = upper_bound
            else:
                upper_bound_particles=[rmax_particles, zmax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_rmin_particles is None) and (bc_zmin_particles is None):
                lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                lower_boundary_conditions_particles = [bc_rmin_particles, bc_zmin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_rmax_particles is None) and (bc_zmax_particles is None):
                upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                upper_boundary_conditions_particles = [bc_rmax_particles, bc_zmax_particles]

        # Sanity check on dimensionality of vector quantities
        assert len(number_of_cells) == 2, Exception('Wrong number of cells specified')
        assert len(lower_bound) == 2, Exception('Wrong number of lower bounds specified')
        assert len(upper_bound) == 2, Exception('Wrong number of upper bounds specified')
        assert len(lower_boundary_conditions) == 2, Exception('Wrong number of lower boundary conditions specified')
        assert len(upper_boundary_conditions) == 2, Exception('Wrong number of upper boundary conditions specified')

        self.number_of_cells = number_of_cells

        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions
        self.lower_bound_particles = lower_bound_particles
        self.upper_bound_particles = upper_bound_particles
        self.lower_boundary_conditions_particles = lower_boundary_conditions_particles
        self.upper_boundary_conditions_particles = upper_boundary_conditions_particles
        self.guard_cells = guard_cells
        self.pml_cells = pml_cells

        self.n_azimuthal_modes = n_azimuthal_modes

        self.moving_window_velocity = moving_window_velocity

        self.refined_regions = refined_regions
        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2,2])
            assert len(region[1]) == 2, Exception('The lo extent of the refined region must be a vector of length 2')
            assert len(region[2]) == 2, Exception('The hi extent of the refined region must be a vector of length 2')
            assert len(region[3]) == 2, Exception('The refinement factor of the refined region must be a vector of length 2')

        self.handle_init(kw)

    def add_refined_region(self, level, lo, hi, refinement_factor=[2,2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_Cartesian2DGrid(_ClassWithInit):
    """
    Two dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, or neumann

    nx: integer
        Number of cells along X (number of nodes=nx+1)

    ny: integer
        Number of cells along Y (number of nodes=ny+1)

    xmin: float
        Position of first node along X [m]

    xmax: float
        Position of last node along X [m]

    ymin: float
        Position of first node along Y [m]

    ymax: float
        Position of last node along Y [m]

    bc_xmin: vector of strings
        Boundary condition at min X: One of periodic, open, dirichlet, or neumann

    bc_xmax: vector of strings
        Boundary condition at max X: One of periodic, open, dirichlet, or neumann

    bc_ymin: vector of strings
        Boundary condition at min Y: One of periodic, open, dirichlet, or neumann

    bc_ymax: vector of strings
        Boundary condition at max Y: One of periodic, open, dirichlet, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 2 specifying the extent of the region,
        and refinement_factor defaulting to [2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    xmin_particles: float, optional
        Position of min particle boundary along X [m]

    xmax_particles: float, optional
        Position of max particle boundary along X [m]

    ymin_particles: float, optional
        Position of min particle boundary along Y [m]

    ymax_particles: float, optional
        Position of max particle boundary along Y [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_xmin_particles: string, optional
        Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal

    bc_xmax_particles: string, optional
        Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal

    bc_ymin_particles: string, optional
        Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal

    bc_ymax_particles: string, optional
        Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx and ny) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 2

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nx=None, ny=None,
                 xmin=None, xmax=None, ymin=None, ymax=None,
                 bc_xmin=None, bc_xmax=None, bc_ymin=None, bc_ymax=None,
                 moving_window_velocity=None, refined_regions=[],lower_bound_particles=None, upper_bound_particles=None,
                 xmin_particles=None, xmax_particles=None, ymin_particles=None, ymax_particles=None,
                 lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None,
                 bc_xmin_particles=None, bc_xmax_particles=None, bc_ymin_particles=None, bc_ymax_particles=None,
                 guard_cells=None, pml_cells=None,
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
            number_of_cells = [nx, ny]
        if lower_bound is None:
            lower_bound = [xmin, ymin]
        if upper_bound is None:
            upper_bound = [xmax, ymax]
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_xmin, bc_ymin]
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_xmax, bc_ymax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (xmin_particles is None) and (ymin_particles is None):
                lower_bound_particles = lower_bound
            else:
                lower_bound_particles = [xmin_particles, ymin_particles]
        if upper_bound_particles is None:
            if (xmax_particles is None) and (ymax_particles is None):
                upper_bound_particles = upper_bound
            else:
                upper_bound_particles=[xmax_particles, ymax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_xmin_particles is None) and (bc_ymin_particles is None):
                lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                lower_boundary_conditions_particles = [bc_xmin_particles, bc_ymin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_xmax_particles is None) and (bc_ymax_particles is None):
                upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                upper_boundary_conditions_particles = [bc_xmax_particles, bc_ymax_particles]

        # Sanity check on dimensionality of vector quantities
        assert len(number_of_cells) == 2, Exception('Wrong number of cells specified')
        assert len(lower_bound) == 2, Exception('Wrong number of lower bounds specified')
        assert len(upper_bound) == 2, Exception('Wrong number of upper bounds specified')
        assert len(lower_boundary_conditions) == 2, Exception('Wrong number of lower boundary conditions specified')
        assert len(upper_boundary_conditions) == 2, Exception('Wrong number of upper boundary conditions specified')
        assert len(lower_bound_particles) == 2, Exception('Wrong number of particle lower bounds specified')
        assert len(upper_bound_particles) == 2, Exception('Wrong number of particle upper bounds specified')
        assert len(lower_boundary_conditions_particles) == 2, Exception('Wrong number of lower particle boundary conditions specified')
        assert len(upper_boundary_conditions_particles) == 2, Exception('Wrong number of upper particle boundary conditions specified')

        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions
        self.lower_bound_particles = lower_bound_particles
        self.upper_bound_particles = upper_bound_particles
        self.lower_boundary_conditions_particles = lower_boundary_conditions_particles
        self.upper_boundary_conditions_particles = upper_boundary_conditions_particles
        self.guard_cells = guard_cells
        self.pml_cells = pml_cells

        self.moving_window_velocity = moving_window_velocity

        self.refined_regions = refined_regions

        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2,2])
            assert len(region[1]) == 2, Exception('The lo extent of the refined region must be a vector of length 2')
            assert len(region[2]) == 2, Exception('The hi extent of the refined region must be a vector of length 2')
            assert len(region[3]) == 2, Exception('The refinement factor of the refined region must be a vector of length 2')

        self.handle_init(kw)

    def add_refined_region(self, level, lo, hi, refinement_factor=[2,2]):
        """Add a refined region.
        level: the refinement level, with 1 being the first level of refinement, 2 being the second etc.
        lo, hi: vectors of length 2 specifying the extent of the region
        refinement_factor: defaulting to [2,2] (relative to next lower level)
        """
        self.refined_regions.append([level, lo, hi, refinement_factor])


class PICMI_Cartesian3DGrid(_ClassWithInit):
    """
    Three dimensional Cartesian grid
    Parameters can be specified either as vectors or separately.
    (If both are specified, the vector is used.)

    Parameters
    ----------
    number_of_cells: vector of integers
        Number of cells along each axis (number of nodes is number_of_cells+1)

    lower_bound: vector of floats
        Position of the node at the lower bound [m]

    upper_bound: vector of floats
        Position of the node at the upper bound [m]

    lower_boundary_conditions: vector of strings
        Conditions at lower boundaries, periodic, open, dirichlet, or neumann

    upper_boundary_conditions: vector of strings
        Conditions at upper boundaries, periodic, open, dirichlet, or neumann

    nx: integer
        Number of cells along X (number of nodes=nx+1)

    ny: integer
        Number of cells along Y (number of nodes=ny+1)

    nz: integer
        Number of cells along Z (number of nodes=nz+1)

    xmin: float
        Position of first node along X [m]

    xmax: float
        Position of last node along X [m]

    ymin: float
        Position of first node along Y [m]

    ymax: float
        Position of last node along Y [m]

    zmin: float
        Position of first node along Z [m]

    zmax: float
        Position of last node along Z [m]

    bc_xmin: string
        Boundary condition at min X: One of periodic, open, dirichlet, or neumann

    bc_xmax: string
        Boundary condition at max X: One of periodic, open, dirichlet, or neumann

    bc_ymin: string
        Boundary condition at min Y: One of periodic, open, dirichlet, or neumann

    bc_ymax: string
        Boundary condition at max Y: One of periodic, open, dirichlet, or neumann

    bc_zmin: string
        Boundary condition at min Z: One of periodic, open, dirichlet, or neumann

    bc_zmax: string
        Boundary condition at max Z: One of periodic, open, dirichlet, or neumann

    moving_window_velocity: vector of floats, optional
        Moving frame velocity [m/s]

    refined_regions: list of lists, optional
        List of refined regions, each element being a list of the format [level, lo, hi, refinement_factor],
        with level being the refinement level, with 1 being the first level of refinement, 2 being the second etc,
        lo and hi being vectors of length 3 specifying the extent of the region,
        and refinement_factor defaulting to [2,2,2] (relative to next lower level)

    lower_bound_particles: vector of floats, optional
        Position of particle lower bound [m]

    upper_bound_particles: vector of floats, optional
        Position of particle upper bound [m]

    xmin_particles: float, optional
        Position of min particle boundary along X [m]

    xmax_particles: float, optional
        Position of max particle boundary along X [m]

    ymin_particles: float, optional
        Position of min particle boundary along Y [m]

    ymax_particles: float, optional
        Position of max particle boundary along Y [m]

    zmin_particles float, optional
        Position of min particle boundary along Z [m]

    zmax_particles: float, optional
        Position of max particle boundary along Z [m]

    lower_boundary_conditions_particles: vector of strings, optional
        Conditions at lower boundaries for particles, periodic, absorbing, reflect or thermal

    upper_boundary_conditions_particles: vector of strings, optional
        Conditions at upper boundaries for particles, periodic, absorbing, reflect or thermal

    bc_xmin_particles: string, optional
        Boundary condition at min X for particles: One of periodic, absorbing, reflect, thermal

    bc_xmax_particles: string, optional
        Boundary condition at max X for particles: One of periodic, absorbing, reflect, thermal

    bc_ymin_particles: string, optional
        Boundary condition at min Y for particles: One of periodic, absorbing, reflect, thermal

    bc_ymax_particles: string, optional
        Boundary condition at max Y for particles: One of periodic, absorbing, reflect, thermal

    bc_zmin_particles: string, optional
        Boundary condition at min Z for particles: One of periodic, absorbing, reflect, thermal

    bc_zmax_particles: string, optional
        Boundary condition at max Z for particles: One of periodic, absorbing, reflect, thermal

    guard_cells: vector of integers, optional
        Number of guard cells used along each direction

    pml_cells: vector of integers, optional
        Number of Perfectly Matched Layer (PML) cells along each direction
    """
    # Note for implementations, as a matter of convenience and flexibility, the user interface allows
    # specifying various quantities using either the individual named attributes (such as nx, ny, and nz) or a
    # vector of values (such as number_of_cells). However, internally, only the vectors are saved and
    # the implementation needs to use the those to access the user input.

    number_of_dimensions = 3

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nx=None, ny=None, nz=None,
                 xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
                 bc_xmin=None, bc_xmax=None, bc_ymin=None, bc_ymax=None, bc_zmin=None, bc_zmax=None,
                 moving_window_velocity=None, refined_regions=[], lower_bound_particles=None, upper_bound_particles=None,
                 xmin_particles=None, xmax_particles=None, ymin_particles=None, ymax_particles=None, zmin_particles=None, zmax_particles=None,
                 lower_boundary_conditions_particles=None, upper_boundary_conditions_particles=None,
                 bc_xmin_particles=None, bc_xmax_particles=None, bc_ymin_particles=None, bc_ymax_particles=None,
                 bc_zmin_particles=None, bc_zmax_particles=None, guard_cells=None, pml_cells=None,
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
            number_of_cells = [nx, ny, nz]
        if lower_bound is None:
            lower_bound = [xmin, ymin, zmin]
        if upper_bound is None:
            upper_bound = [xmax, ymax, zmax]
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_xmin, bc_ymin, bc_zmin]
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_xmax, bc_ymax, bc_zmax]

        # Sanity check and init of input arguments related to particle boundary parameters
        # By default, if not specified, particle boundary values are the same as field boundary values
        # By default, if not specified, particle boundary conditions are the same as field boundary conditions
        if lower_bound_particles is None:
            if (xmin_particles is None) and (ymin_particles is None) and (zmin_particles is None):
                lower_bound_particles = lower_bound
            else:
                lower_bound_particles = [xmin_particles, ymin_particles, zmin_particles]
        if upper_bound_particles is None:
            if (xmax_particles is None) and (ymax_particles is None) and (zmax_particles is None):
                upper_bound_particles = upper_bound
            else:
                upper_bound_particles = [xmax_particles, ymax_particles, zmax_particles]

        if lower_boundary_conditions_particles is None:
            if (bc_xmin_particles is None) and (bc_ymin_particles is None) and (bc_zmin_particles is None):
                lower_boundary_conditions_particles = lower_boundary_conditions
            else:
                lower_boundary_conditions_particles = [bc_xmin_particles, bc_ymin_particles, bc_zmin_particles]
        if upper_boundary_conditions_particles is None:
            if (bc_xmax_particles is None) and (bc_ymax_particles is None) and (bc_zmax_particles is None):
                upper_boundary_conditions_particles = upper_boundary_conditions
            else:
                upper_boundary_conditions_particles = [bc_xmax_particles, bc_ymax_particles, bc_zmax_particles]

        # Sanity check on number of arguments of vector quantities
        assert len(number_of_cells) == 3, Exception('Wrong number of cells specified')
        assert len(lower_bound) == 3, Exception('Wrong number of lower bounds specified')
        assert len(upper_bound) == 3, Exception('Wrong number of upper bounds specified')
        assert len(lower_boundary_conditions) == 3, Exception('Wrong number of lower boundary conditions specified')
        assert len(upper_boundary_conditions) == 3, Exception('Wrong number of upper boundary conditions specified')
        assert len(lower_bound_particles) == 3, Exception('Wrong number of particle lower bounds specified')
        assert len(upper_bound_particles) == 3, Exception('Wrong number of particle upper bounds specified')
        assert len(lower_boundary_conditions_particles) == 3, Exception('Wrong number of particle lower boundary conditions specified')
        assert len(upper_boundary_conditions_particles) == 3, Exception('Wrong number of particle upper boundary conditions specified')

        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions
        self.lower_bound_particles = lower_bound_particles
        self.upper_bound_particles = upper_bound_particles
        self.lower_boundary_conditions_particles = lower_boundary_conditions_particles
        self.upper_boundary_conditions_particles = upper_boundary_conditions_particles
        self.guard_cells = guard_cells
        self.pml_cells = pml_cells

        self.moving_window_velocity = moving_window_velocity

        self.refined_regions = refined_regions
        for region in self.refined_regions:
            if len(region) == 3:
                region.append([2,2,2])
            assert len(region[1]) == 3, Exception('The lo extent of the refined region must be a vector of length 3')
            assert len(region[2]) == 3, Exception('The hi extent of the refined region must be a vector of length 3')
            assert len(region[3]) == 3, Exception('The refinement factor of the refined region must be a vector of length 3')

        self.handle_init(kw)

    def add_refined_region(self, level, lo, hi, refinement_factor=[2,2,2]):
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
