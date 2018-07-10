"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys


# Physics objects
# ---------------


class PICMI_Species(object):
    """
    Species
      - type: a string specifying an elementary particle, atom, or other, as defined in the openPMD species type extension
      - name: name of the species (optional)
      - charge_state: charge state of the species (applies to atoms) (optional)
      - charge: particle charge, required when type is not specified, otherwise determined from type [C]
      - mass: particle mass, required when type is not specified, otherwise determined from type [kg]
    """

    def __init__(self, type=None, name=None, charge_state=None, charge=None, mass=None,
                 **kw):

        self.type = type
        self.name = name
        self.charge = charge
        self.charge_state = charge_state
        self.mass = mass

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


class PICMI_GaussianBeam(object):
    """
    Describes a Gaussian distribution of particles
      - species: Particle species
      - number_real_particles: Number of real particles in the beam.
      - number_sim_particles: Number of simulation particles in the beam.
      - T0=0.: Time at which parameters are specified [s]
      - Xmean=0.: Mean X position [m]
      - Ymean=0.: Mean Y position [m]
      - Zmean=0.: Mean Z position [m]
      - Xrms=0.: R.M.S. size along X [m]
      - Yrms=0.: R.M.S. size along Y [m]
      - Zrms=0.: R.M.S. size along Z [m]
      - UXmean=0.: Mean velocity (gamma*V) along X [m/s]
      - UYmean=0.: Mean velocity (gamma*V) along X [m/s]
      - UZmean=0.: Mean velocity (gamma*V) along X [m/s]
      - UXrms=0.: R.M.S. velocity (gamma*V) spread along X [m/s]
      - UYrms=0.: R.M.S. velocity (gamma*V) spread along Y [m/s]
      - UZrms=0.: R.M.S. velocity (gamma*V) spread along Z [m/s]
      - UXdiv=0.: Velocity (gamma*V) divergence along X [m/s/m]
      - UYdiv=0.: Velocity (gamma*V) divergence along Y [m/s/m]
      - UZdiv=0.: Velocity (gamma*V) divergence along Z [m/s/m]
      - density_func=None: Function modulating density as a function of x, y, z and/or time
      - array_func=None: Array modulating density as a function of x, y, z and/or time
    """
    def __init(self, species, number_real_particles, number_sim_particles, T0=0.,
               Xmean=0., Ymean=0., Zmean=0.,
               Xrms=0., Yrms=0., Zrms=0.,
               UXmean=0., UYmean=0., UZmean=0.,
               UXrms=0., UYrms=0., UZrms=0.,
               UXdiv=0., UYdiv=0., UZdiv=0.,
               density_func=None, array_func=None,
               **kw):
        self.species = species
        self.number_real_particles = number_real_particles
        self.number_sim_particles = number_sim_particles
        self.T0 = T0
        self.Xmean = Xmean
        self.Ymean = Ymean
        self.Zmean = Zmean
        self.Xrms = Xrms
        self.Yrms = Yrms
        self.Zrms = Zrms
        self.UXmean = UXmean
        self.UYmean = UYmean
        self.UZmean = UZmean
        self.UXrms = UXrms
        self.UYrms = UYrms
        self.UZrms = UZrms
        self.density_func = density_func
        self.array_func = array_func

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


class PICMI_Plasma(object):
    """
    Describes a uniform density plasma
      - species: Particle species or list of species
      - density: Plasma density [m^-3].
      - xmin=-infinity: Min position of box along X [m]
      - xmax=+infinity: Max position of box along X [m]
      - ymin=-infinity: Min position of box along Y [m]
      - ymax=+infinity: Max position of box along Y [m]
      - zmin=-infinity: Min position of box along Z [m]
      - zmax=+infinity: Max position of box along Z [m]
      - vthx=0.: Thermal velocity along X [m/s]
      - vthy=0.: Thermal velocity along Y [m/s]
      - vthz=0.: Thermal velocity along Z [m/s]
      - vxmean=0.: Mean velocity along X [m/s]
      - vymean=0.: Mean velocity along Y [m/s]
      - vzmean=0.: Mean velocity along Z [m/s]
      - number_per_cell: Number of particles per cell (randomly placed)
                         Only one of number_per_cell or number_per_cell_each_dim should be specified.
      - number_per_cell_each_dim: Number of particles along each axis (for regularly placed particles)
                                  Only one of number_per_cell or number_per_cell_each_dim should be specified.
      - density_func=None: Function modulating density as a function of x, y, z and/or time.
      - array_func=None: Array modulating density as a function of x, y, z and/or time.
      - fill_in=False: Flags whether to fill in the empty spaced opened up when the grid moves

    """

    def __init__(self, species, density,
                 xmin=None, xmax=None,
                 ymin=None, ymax=None,
                 zmin=None, zmax=None,
                 vthx=0., vthy=0., vthz=0., vxmean=0., vymean=0., vzmean=0.,
                 number_per_cell=None, number_per_cell_each_dim=None, density_func=None, array_func=None,
                 fill_in=False,
                 **kw):
        if not isinstance(species, list):
            species = [species]
        self.species = species
        self.density = density
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.vthx = vthx
        self.vthy = vthy
        self.vthz = vthz
        self.vxmean = vxmean
        self.vymean = vymean
        self.vzmean = vzmean
        self.number_per_cell = number_per_cell
        self.number_per_cell_each_dim = number_per_cell_each_dim
        self.density_func = density_func
        self.array_func = array_func
        self.fill_in = fill_in

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


class PICMI_ParticleList(object):
    """
    Load particles at the specified positions and velocities
      - species: Particle species
      - weight: Particle weight, number of real particles per simulation particle
      - x=0.: List of x positions of the particles [m]
      - y=0.: List of y positions of the particles [m]
      - z=0.: List of z positions of the particles [m]
      - ux=0.: List of ux positions of the particles (ux = gamma*vx) [m/s]
      - uy=0.: List of uy positions of the particles (uy = gamma*vy) [m/s]
      - uz=0.: List of uz positions of the particles (uz = gamma*vz) [m/s]
    """
    def __init__(self, x=0., y=0., z=0., ux=0., uy=0., uz=0.,
                 **kw):
        # --- Get length of arrays, set to one for scalars
        lenx = np.size(x)
        leny = np.size(y)
        lenz = np.size(z)
        lenux = np.size(ux)
        lenuy = np.size(uy)
        lenuz = np.size(uz)

        maxlen = max(lenx, leny, lenz, lenux, lenuy, lenuz)
        assert lenx==maxlen or lenx==1, "Length of x doesn't match len of others"
        assert leny==maxlen or leny==1, "Length of y doesn't match len of others"
        assert lenz==maxlen or lenz==1, "Length of z doesn't match len of others"
        assert lenux==maxlen or lenux==1, "Length of ux doesn't match len of others"
        assert lenuy==maxlen or lenuy==1, "Length of uy doesn't match len of others"
        assert lenuz==maxlen or lenuz==1, "Length of uz doesn't match len of others"

        if lenx == 1:
            x = np.array(x)*np.ones(maxlen)
        if leny == 1:
            y = np.array(y)*np.ones(maxlen)
        if lenz == 1:
            z = np.array(z)*np.ones(maxlen)
        if lenux == 1:
            ux = np.array(ux)*np.ones(maxlen)
        if lenuy == 1:
            uy = np.array(uy)*np.ones(maxlen)
        if lenuz == 1:
            uz = np.array(uz)*np.ones(maxlen,'d')

        self.species = species
        self.weight = weight
        self.x = x
        self.y = y
        self.z = z
        self.ux = ux
        self.uy = uy
        self.uz = uz

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


# ------------------
# Numerics Objects
# ------------------


class PICMI_ParticleDistributionInjector(object):
    """
    Describes the injection of particles from a plane
      - species: Particle species to inject
      - distribution: Particle distribution to inject
      - method: InPlace - method of injection. One of 'InPlace', or 'Plane'
      - X0=0.: Position of the particle centroid in X [m]
      - Y0=0.: Position of the particle centroid in Y [m]
      - Z0=0.: Position of the particle centroid in Z [m]
      - Xplane=0.: Position of the plane of injection in X [m]
      - Yplane=0.: Position of the plane of injection in Y [m]
      - Zplane=0.: Position of the plane of injection in Z [m]
      - VXplane=0.: Velocity of the plane of injection in X [m/s]
      - VYplane=0.: Velocity of the plane of injection in Y [m/s]
      - VZplane=0.: Velocity of the plane of injection in Z [m/s]
      - XVecPlane=0.: Component along X of vector normal to injection plane
      - YVecPlane=0.: Component along Y of vector normal to injection plane
      - ZVecPlane=1.: Component along Z of vector normal to injection plane
    """
    def __init__(self, species, distribution, method='InPlace', X0=0., Y0=0., Z0=0.,
                 Xplane=0., Yplane=0., Zplane=0.,
                 VXplane=0., VYplane=0., VZplane=0.,
                 XVecPlane=0., YVecPlane=0., ZVecPlane=1.,
                 **kw):
        self.species = species
        self.distribution = distribution
        self.method = method
        self.X0 = X0
        self.Y0 = Y0
        self.Z0 = Z0
        self.Xplane = Xplane
        self.Yplane = Yplane
        self.Zplane = Zplane
        self.VXplane = VXplane
        self.VYplane = VYplane
        self.VZplane = VZplane
        self.XVecPlane = XVecPlane
        self.YVecPlane = YVecPlane
        self.ZVecPlane = ZVecPlane

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


class PICMI_Grid(object):
    """
    Grid
      - nx: Number of cells along X (Nb nodes=nx+1)
      - ny: Number of cells along Y (Nb nodes=ny+1)
      - nr: Number of cells along R (Nb nodes=nr+1)
      - nz: Number of cells along Z (Nb nodes=nz+1)
      - nm: Number of azimuthal modes
      - xmin: Position of first node along X [m]
      - xmax: Position of last node along X [m]
      - ymin: Position of first node along Y [m]
      - ymax: Position of last node along Y [m]
      - rmax: Position of last node along R [m]
      - zmin: Position of first node along Z [m]
      - zmax: Position of last node along Z [m]
      - bcxmin: Boundary condition at min X: One of periodic, open, dirichlet, or neumann
      - bcxmax: Boundary condition at max X: One of periodic, open, dirichlet, or neumann
      - bcymin: Boundary condition at min Y: One of periodic, open, dirichlet, or neumann
      - bcymax: Boundary condition at max Y: One of periodic, open, dirichlet, or neumann
      - bcrmax: Boundary condition at max R: One of open, dirichlet, or neumann
      - bczmin: Boundary condition at min Z: One of periodic, open, dirichlet, or neumann
      - bczmax: Boundary condition at max Z: One of periodic, open, dirichlet, or neumann
      - moving_window_velocity: An array of the moving frame velocity in each direction [m/s]
    """

    def __init__(self, nx=None, ny=None, nr=None, nz=None, nm=None,
                 xmin=None, xmax=None, ymin=None, ymax=None, rmax=None, zmin=None, zmax=None,
                 bcxmin=None, bcxmax=None, bcymin=None, bcymax=None, bcrmax=None, bczmin=None, bczmax=None,
                 moving_window_velocity=None,  **kw):
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
        self.moving_window_velocity = moving_window_velocity

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError

    def getmins(self, **kw):
        raise NotImplementedError

    def getmaxs(self, **kw):
        raise NotImplementedError

    def getdims(self, **kw):
        raise NotImplementedError


class PICMI_EM_solver(object):
    """
    EM_solver
      - method: One of Yee, CK, CKC, Lehe, PSTD, PSATD, or GPSTD
      - norderx: Order of stencil in X (-1=infinite)
      - nordery: Order of stencil in Y (-1=infinite)
      - norderr: Order of stencil in R (-1=infinite)
      - norderz: Order of stencil in Z (-1=infinite)
      - l_nodal: Quantities are at nodes if True, staggered otherwise

    Methods:
      - add: Add object to solver (e.g. laser)
    """

    methods_list = ['Yee', 'CK', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD']

    def __init__(self, method=None,
                 norderx=None, nordery=None, norderr=None, norderz=None,
                 l_nodal=None, **kw):

        assert method is None or method in EM_solver.Methods_list, Exception('method has incorrect value')

        self.method = method
        self.norderx = norderx
        self.nordery = nordery
        self.norderr = norderr
        self.norderz = norderz
        self.l_nodal = l_nodal

        self.objects = []

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError

    def add(self, object):
        self.objects.append(object)


class PICMI_ES_solver(object):
    """
    ES_solver
      - method: One of FFT, or Multigrid
    """

    methods_list = ['FFT', 'Multigrid']

    def __init__(self, method=None):

        assert method is None or method in ES_solver.methods_list, Exception('method has incorrect value')

        self.method = method

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


class PICMI_Simulation(object):
    """
    Simulation
      - timestep=0.: Absolute time step size of the simulation [s]
                     (set to 0 to have the timestep relative to the CFL limit)
      - timestep_over_cfl=1.: Ratio of the time step size to the Courant-Friedrich-Lewy limit
                              (used only if timestep is 0; should raise an error when the code does not have a well-defined CFL)
      - max_step: Maximum number of time steps
      - max_time: Maximum time to run the simulation [s]
      - verbose: Verbosity flag
    """

    def __init__(self, timestep=0., timestep_over_cfl=1., max_step=None, max_time=None, verbose=None,
                 **kw):
        self.timestep = timestep
        self.timestep_over_cfl = timestep_over_cfl
        self.verbose = verbose
        self.max_step = max_step
        self.max_time = max_time

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError

    def step(self, nsteps=1):
        raise NotImplementedError


# -----------------------------------------
# Injecting a laser pulse in the Simulation
# -----------------------------------------


def add_laser_pulse(simulation, laser_profile, injection_method):
    """
    Specify laser pulses that are injected in the simulation
      - simulation: a Simulation object, which contains all the relevant data for a simulation.
      - laser_profile: one of laser profile objects
                       Specifies the **physical** properties of the laser pulse.
                       (e.g. spatial and temporal profile, wavelength, amplitude, etc.)
      - injection_method: a laser injector object (optional)
                          Specifies how the laser is injected (numerically) into the simulation
                          (e.g. through a laser antenna, or directly added to the mesh).
                          This argument describes an **algorithm**, not a physical object.
                          It is optional. (It is up to each code to define the default method
                          of injection, if the user does not provide injection_method)
    """
    raise NotImplementedError


# Laser profile objects


class PICMI_Gaussian_laser(object):
    """
    Gaussian Laser
      - wavelength: Laser wavelength.
      - waist: Waist of the Gaussian pulse at focus [m].
      - duration: Duration of the Gaussian pulse [s].
      - focal_position: Position of the laser focus. [m]
      - x0: Position of the laser centroid in X at time 0. [m]
      - y0: Position of the laser centroid in Y at time 0. [m]
      - z0: Position of the laser centroid in Z at time 0. [m]
      - pol_angle: Angle of polarization (relative to X). [radians]
      - a0: Normalized vector potential at focus.
            Specify either a0 or E0 (E0 takes precedence)
      - E0: Maximum amplitude of the laser field. [V/m]
            Specify either a0 or E0 (E0 takes precedence)
    """
    def __init__(self, wavelength, waist, duration, focal_position=0., x0=0., y0=0., z0=0.,
                 pol_angle=None, a0=None, E0=None,
                 **kw):
        k0 = 2.*math.pi/wavelength
        if E0 is None:
            from scipy import constants
            E0 = a0*constants.electron_mass*constants.speed_of_light**2*k0/constants.elementary_charge
        if a0 is None:
            from scipy import constants
            a0 = E0/(constants.electron_mass*constants.speed_of_light**2*k0/constants.elementary_charge)

        self.wavelength = wavelength
        self.k0 = k0
        self.waist = waist
        self.duration = duration
        self.focal_position = focal_position
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.pol_angle = pol_angle
        self.a0 = a0
        self.E0 = E0

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError


# Laser injection methods


class PICMI_Laser_antenna(object):
    """
    Laser antenna injection method
      - laser: Laser object to be injected
      - antenna_x0=0.: Position of antenna launching the laser along X. [m]
      - antenna_y0=0.: Position of antenna launching the laser along Y. [m]
      - antenna_z0=0.: Position of antenna launching the laser along Z. [m]
      - antenna_xvec=0.: Component along X of vector normal to antenna plane.
      - antenna_yvec=0.: Component along Y of vector normal to antenna plane.
      - antenna_zvec=1.: Component along Z of vector normal to antenna plane.
    """
    def __init__(self, laser, antenna_x0=0., antenna_y0=0., antenna_z0=0.,
                 antenna_xvec=0., antenna_yvec=0., antenna_zvec=1., **kw):

        self.laser = laser
        self.antenna_x0 = antenna_x0
        self.antenna_y0 = antenna_y0
        self.antenna_z0 = antenna_z0
        self.antenna_xvec = antenna_xvec
        self.antenna_yvec = antenna_yvec
        self.antenna_zvec = antenna_zvec

        self.init(**kw)

    def init(self, **kw):
        raise NotImplementedError
