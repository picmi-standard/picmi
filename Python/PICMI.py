"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
from CODATA2014 import *


class PICMI_Grid(object):
    """
    Grid
      - nx: Number of cells along X (Nb nodes=nx+1)
      - ny: Number of cells along Y (Nb nodes=ny+1)
      - nr: Number of cells along R (Nb nodes=nr+1)
      - nz: Number of cells along Z (Nb nodes=nz+1)
      - nm: Number of azimuthal modes
      - xmin: Position of first node along X
      - xmax: Position of last node along X
      - ymin: Position of first node along Y
      - ymax: Position of last node along Y
      - rmax: Position of last node along R
      - zmin: Position of first node along Z
      - zmax: Position of last node along Z
      - bcxmin: Boundary condition at min X: periodic/open/dirichlet/neumann
      - bcxmax: Boundary condition at max X: periodic/open/dirichlet/neumann
      - bcymin: Boundary condition at min Y: periodic/open/dirichlet/neumann
      - bcymax: Boundary condition at max Y: periodic/open/dirichlet/neumann
      - bcrmax: Boundary condition at max R: open/dirichlet/neumann
      - bczmin: Boundary condition at min Z: periodic/open/dirichlet/neumann
      - bczmax: Boundary condition at max Z: periodic/open/dirichlet/neumann
    """

    def __init__(self, nx=None, ny=None, nr=None, nz=None, nm=None,
                 xmin=None, xmax=None, ymin=None, ymax=None, rmax=None, zmin=None, zmax=None,
                 bcxmin=None, bcxmax=None, bcymin=None, bcymax=None, bcrmax=None, bczmin=None, bczmax=None,
                 **kw):
        self.nx = nx
        self.ny = ny
        self.nr = nr
        self.nz = nz
        self.nm = nm
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.rmax = rmax
        self.zmin = zmin
        self.zmax = zmax
        self.bcxmin = bcxmin
        self.bcxmax = bcxmax
        self.bcymin = bcymin
        self.bcymax = bcymax
        self.bcrmax = bcrmax
        self.bczmin = bczmin
        self.bczmax = bczmax

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


class PICMI_EM_solver(object):
    """
    EM_solver
      - method: Yee/CK/CKC/Lehe/PSTD/PSATD/GPSTD
      - norderx: Order of stencil in X (-1=infinite)
      - nordery: Order of stencil in Y (-1=infinite)
      - norderr: Order of stencil in R (-1=infinite)
      - norderz: Order of stencil in Z (-1=infinite)
      - l_nodal: Quantities are at nodes if True, staggered otherwise
    """

    Methods_list = ['Yee', 'CK', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD']

    def __init__(self, Method=None,
                 norderx=None, nordery=None, norderr=None, norderz=None,
                 l_nodal=None,
                 current_deposition_algo=None, charge_deposition_algo=None,
                 field_gathering_algo=None, particle_pusher_algo=None, **kw):

        assert Method is None or Method in EM_solver.Methods_list, Exception('Method has incorrect value')

        self.norderx = norderx
        self.nordery = nordery
        self.norderr = norderr
        self.norderz = norderz
        self.l_nodal = l_nodal
        self.current_deposition_algo = current_deposition_algo
        self.charge_deposition_algo = charge_deposition_algo
        self.field_gathering_algo = field_gathering_algo
        self.particle_pusher_algo = particle_pusher_algo

        self.init(**kw)
    
    def init(self, **kw):
        raise NotImplementedError


class PICMI_Species(object):
    """
    Species
      - type: an elementary particle or atom, or user-defined type
      - name: name of the species
      - sid: unique identification number of the species
      - charg_state: charge state of the species (applies to atoms)
      - charge: (multiplied by Charge_state if atom)
      - mass: mass of species particle
      - weight: weight of the species
    """

    def __init__(self, type=None, name=None, sid=None,
                 charge_state=None, charge=None, mass=None,
                 weight=1., **kw):

        self.type = type
        self.name = name
        self.sid = sid
        self.charge_state = charge_state
        self.charge = charge
        self.mass = mass
        self.weight = weight

        self.init(**kw)

    def PICMI_init_charge(self, charge, charge_state):
        """
        The input argument, charge, takes precedent, followed by charge_state, then
        the possible charge of the particle type, and finally defaulting to zero.
        """
        self.charge_state = charge_state
        if charge is not None:
            self.charge = charge
        elif charge_state is not None:
            self.charge = echarge*charge_state
        else:
            try:
                self.charge = self.type.charge
            except AttributeError:
                self.charge = 0.

    def PICMI_init_mass(self, mass):
        """
        The input argument, mass, takes precedent, followed by the mass of the particle
        and finally defaulting to zero.
        """
        if mass is not None:
            self.mass = mass
        else:
            try:
                self.mass = self.type.mass
            except AttributeError:
                self.mass = 0.

    def init(self, **kw):
        raise NotImplementedError

    def add_particles(self, n=None,
                      x=None, y=None, z=None,
                      ux=None, uy=None, uz=None, w=None,
                      unique_particles=None, **kw):
        raise NotImplementedError


class PICMI_Simulation(object):
    """
    Simulation
      - plot_in: Diagnostic output interval
      - verbose: Verbosity flag
      - cfl: Courant-Friedrich-Lewy limit
    """

    def __init__(self, plot_int=None, verbose=None, cfl=None, **kw):
        self.plot_int = plot_int
        self.verbose = verbose
        self.cfl = cfl

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError

    def step(self, nsteps=1):
        raise NotImplementedError
        
    def finalize(self):
        raise NotImplementedError
