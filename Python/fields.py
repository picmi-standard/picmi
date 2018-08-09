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
      - grid: grid object to be used by the solver
      - method: One of Yee, CKC, Lehe, PSTD, PSATD, or GPSTD
      - stencil_order: Order of stencil for each axis (-1=infinite) (vector)
      - cfl=None: Fraction of the Courant-Friedrich-Lewy criteria [1]
      - l_nodal=None: Quantities are at nodes if True, staggered otherwise
      - source_smoother=None: Smoother to apply to the sources
      - field_smoother=None: Smoother to apply to the fields
    """

    methods_list = ['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD']

    def __init__(self, grid, method=None, stencil_order=None, cfl=None, l_nodal=None,
                 source_smoother=None, field_smoother=None,
                 **kw):

        assert method is None or method in PICMI_ElectromagneticSolver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectromagneticSolver.methods_list))

        self.grid = grid
        self.method = method
        self.cfl = cfl
        self.stencil_order = stencil_order
        self.l_nodal = l_nodal
        self.source_smoother = source_smoother
        self.field_smoother = field_smoother

        self.handle_init(kw)


class PICMI_Electrostatic_solver(_ClassWithInit):
    """
    Electrostatic field solver
      - grid: grid object to be used by the solver
      - method: One of FFT, or Multigrid
    """

    methods_list = ['FFT', 'Multigrid']

    def __init__(self, grid, method=None):

        assert method is None or method in PICMI_Electrostatic_solver.methods_list, \
               Exception('method must be one of '+', '.join(PICMI_Electrostatic_solver.methods_list))

        self.grid = grid
        self.method = method

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


class PICMI_BinomialSmoother(_ClassWithInit):
    """
    Descibes a binomial smoother operator (applied to grids)
    - n_pass: Number of passes along each axis (vector)
    - compensation: Flags whether to apply comensation (logical)
    """
    def __init__(self, n_pass=None, compensation=None, **kw):
        self.n_pass = n_pass
        self.compensation = compensation

        self.handle_init(kw)


class PICMI_CylindricalGrid(_ClassWithInit):
    """
    Axisymmetric, cylindrical grid
    Parameters can be specified either as vectors or separately.

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

      - moving_window_velocity: Moving frame Z velocity [m/s]
    """

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nr=None, nz=None, n_azimuthal_modes=None,
                 rmin=None, rmax=None, zmin=None, zmax=None,
                 bc_rmin=None, bc_rmax=None, bc_zmin=None, bc_zmax=None,
                 moving_window_velocity=None,  **kw):

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
            number_of_cells = [nr, nz]
        else:
            nr, nz = number_of_cells
        if lower_bound is None:
            lower_bound = [rmin, zmin]
        else:
            rmin, zmin = lower_bound
        if upper_bound is None:
            upper_bound = [rmax, zmax]
        else:
            rmax, zmax = upper_bound
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_rmin, bc_zmin]
        else:
            bc_rmin, bc_zmin = lower_boundary_conditions
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_rmax, bc_zmax]
        else:
            bc_rmax, bc_zmax = upper_boundary_conditions

        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions

        self.nr = nr
        self.nz = nz
        self.n_azimuthal_modes = n_azimuthal_modes
        self.rmin = rmin
        self.rmax = rmax
        self.zmin = zmin
        self.zmax = zmax
        self.bc_rmin = bc_rmin
        self.bc_rmax = bc_rmax
        self.bc_zmin = bc_zmin
        self.bc_zmax = bc_zmax

        self.moving_window_velocity = moving_window_velocity

        self.handle_init(kw)


class PICMI_Cartesian2DGrid(_ClassWithInit):
    """
    Two dimensional Carteisan grid
    Parameters can be specified either as vectors or separately.

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
    """

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nx=None, ny=None,
                 xmin=None, xmax=None, ymin=None, ymax=None,
                 bc_xmin=None, bc_xmax=None, bc_ymin=None, bc_ymax=None,
                 moving_window_velocity=None,
                 **kw):

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
        else:
            nx, ny = number_of_cells
        if lower_bound is None:
            lower_bound = [xmin, ymin]
        else:
            xmin, ymin = lower_bound
        if upper_bound is None:
            upper_bound = [xmax, ymax]
        else:
            xmax, ymax = upper_bound
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_xmin, bc_ymin]
        else:
            bc_xmin, bc_ymin = lower_boundary_conditions
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_xmax, bc_ymax]
        else:
            bc_xmax, bc_ymax = upper_boundary_conditions

        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions

        self.nx = nx
        self.ny = ny
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.bc_xmin = bc_xmin
        self.bc_xmax = bc_xmax
        self.bc_ymin = bc_ymin
        self.bc_ymax = bc_ymax

        self.moving_window_velocity = moving_window_velocity

        self.handle_init(kw)


class PICMI_Cartesian3DGrid(_ClassWithInit):
    """
    Three dimensional Carteisan grid
    Parameters can be specified either as vectors or separately.

      - number_of_cells: Number of cells along each axis (number of nodes is number_of_cells+1) (vector)
      - lower_bound: Position of the node at the lower bound (vector) [m]
      - upper_bound: Position of the node at the upper bound (vector) [m]
      - lower_boundary_conditions: Conditions at lower boundaries, periodic, open, dirichlet, or neumann (vector)
      - upper_boundary_conditions: Conditions at upper boundaries, periodic, open, dirichlet, or neumann (vector)

      - nx: Number of cells along X (number of nodes=nx+1)
      - ny: Number of cells along Y (number of nodes=ny+1)
      - nz: Number of cells along Z (number of nodes=nz+1)
      - xmin: Position of first node along X [m]
      - xmax: Position of last node along X [m]
      - ymin: Position of first node along Y [m]
      - ymax: Position of last node along Y [m]
      - zmin: Position of first node along Z [m]
      - zmax: Position of last node along Z [m]
      - bc_xmin: Boundary condition at min X: One of periodic, open, dirichlet, or neumann
      - bc_xmax: Boundary condition at max X: One of periodic, open, dirichlet, or neumann
      - bc_ymin: Boundary condition at min Y: One of periodic, open, dirichlet, or neumann
      - bc_ymax: Boundary condition at max Y: One of periodic, open, dirichlet, or neumann
      - bc_zmin: Boundary condition at min Z: One of periodic, open, dirichlet, or neumann
      - bc_zmax: Boundary condition at max Z: One of periodic, open, dirichlet, or neumann

      - moving_window_velocity: Moving frame velocity (vector) [m/s]
    """

    def __init__(self, number_of_cells=None, lower_bound=None, upper_bound=None,
                 lower_boundary_conditions=None, upper_boundary_conditions=None,
                 nx=None, ny=None, nz=None,
                 xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
                 bc_xmin=None, bc_xmax=None, bc_ymin=None, bc_ymax=None, bc_zmin=None, bc_zmax=None,
                 moving_window_velocity=None,
                 **kw):

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
        else:
            nx, ny, nz = number_of_cells
        if lower_bound is None:
            lower_bound = [xmin, ymin, zmin]
        else:
            xmin, ymin, zmin = lower_bound
        if upper_bound is None:
            upper_bound = [xmax, ymax, zmax]
        else:
            xmax, ymax, zmax = upper_bound
        if lower_boundary_conditions is None:
            lower_boundary_conditions = [bc_xmin, bc_ymin, bc_zmin]
        else:
            bc_xmin, bc_ymin, bc_zmin = lower_boundary_conditions
        if upper_boundary_conditions is None:
            upper_boundary_conditions = [bc_xmax, bc_ymax, bc_zmax]
        else:
            bc_xmax, bc_ymax, bc_zmax = upper_boundary_conditions

        self.number_of_cells = number_of_cells
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.lower_boundary_conditions = lower_boundary_conditions
        self.upper_boundary_conditions = upper_boundary_conditions

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.bc_xmin = bc_xmin
        self.bc_xmax = bc_xmax
        self.bc_ymin = bc_ymin
        self.bc_ymax = bc_ymax
        self.bc_zmin = bc_zmin
        self.bc_zmax = bc_zmax

        self.moving_window_velocity = moving_window_velocity

        self.handle_init(kw)
