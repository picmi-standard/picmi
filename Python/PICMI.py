"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys

code = None

class _ClassWithInit(object):
    def init(self, **kw):
        pass


# Physics objects
# ---------------


class PICMI_Species(_ClassWithInit):
    """
    Species
      - particle_type=None: A string specifying an elementary particle, atom, or other, as defined in the openPMD 2 species type extension
      - name=None: Name of the species
      - charge_state=None: Charge state of the species (applies to atoms) [1]
      - charge=None: Particle charge, required when type is not specified, otherwise determined from type [C]
      - mass=None: Particle mass, required when type is not specified, otherwise determined from type [kg]
      - initial_distribution=None: The initial distribution loaded at t=0. Must be one of the standard distributions implemented.
      - particle_shape: Particle shape used for deposition and gather ; if None, the default from the `Simulation` object will be used. Possible values are 'NGP', 'linear', 'quadratic', 'cubic'
    """

    def __init__(self, particle_type=None, name=None, charge_state=None, charge=None, mass=None,
                 initial_distribution=None, particle_shape=None, **kw):
        self.particle_type = particle_type
        self.name = name
        self.charge = charge
        self.charge_state = charge_state
        self.mass = mass
        self.initial_distribution = initial_distribution

        self.interactions = []

        self.init(**kw)

    def activate_ionization(self, model, target_species):
        # --- TODO: One way of handling interactions is to add a class for each type
        # ---       of interaction. Instances would be added to the interactions list
        # ---       instead of the list of parameters.
        self.interactions.append(['ionization', model, target_species])


class PICMI_MultiSpecies(_ClassWithInit):
    """
    INCOMPLETE: proportions argument is not implemented
    Multiple species that are initialized with the same distribution
      - particle_types: A list of strings specifying for each species an elementary particle, atom, or other, as defined in the openPMD 2 species type extension
      - names: Names of the species (optional)
      - charge_states: Charge states of the species (applies to atoms) (optional)
      - charges: Particle charges, required when type is not specified, otherwise determined from type [C]
      - masses: Particle masses, required when type is not specified, otherwise determined from type [kg]
      - proportions: Proportions of the initial distribution made up by each species
    """
    # --- Note to developer: This class attribute needs to be set to the Species class
    # --- defined in the codes PICMI implementation.
    Species_class = None

    def __init__(self, particle_types=[], names=[], charge_states=[], charges=[], masses=[],
                 initial_distribution=None, proportions=[],
                 **kw):

        self.particle_types = particle_types
        self.names = names
        self.charges = charges
        self.charge_states = charge_states
        self.masses = masses
        self.initial_distribution = initial_distribution
        self.proportions = proportions

        self.nspecies = 0
        self.check_nspecies(particle_types)
        self.check_nspecies(names)
        self.check_nspecies(charges)
        self.check_nspecies(charge_states)
        self.check_nspecies(masses)
        self.check_nspecies(proportions)

        # --- Create the instances of each species
        self.species_instances_list = []
        self.species_instances_dict = {}
        for i in range(self.nspecies):
            if particle_types:
                particle_type = particle_types[i]
            else:
                particle_type = None
            if names:
                name = names[i]
            else:
                name = None
            if charges:
                charge = charges[i]
            else:
                charge = None
            if charge_states:
                charge_state = charge_states[i]
            else:
                charge_state = None
            if masses:
                mass = masses[i]
            else:
                mass = None
            if proportions:
                proportion = proportions[i]
            else:
                proportion = None
            specie = PICMI_MultiSpecies.Species_class(particle_type = particle_type,
                                                      name = name,
                                                      charge = charge,
                                                      charge_state = charge_state,
                                                      mass = mass,
                                                      initial_distribution = initial_distribution)
                                                      #proportion = proportion)
            self.species_instances_list.append(specie)
            if name is not None:
                self.species_instances_dict[name] = specie

        self.init(**kw)

    def check_nspecies(self, var):
        if len(var) > 0:
            assert self.nspecies == 0 or self.nspecies == len(var), Exception('All inputs must have the same length')
            self.nspecies = len(var)

    def __len__(self):
        return self.nspecies

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.species_instances_dict[key]
        else:
            return self.species_instances_list[key]


class PICMI_GaussianBunchDistribution(_ClassWithInit):
    """
    Describes a Gaussian distribution of particles
      - n_physical_particles: Number of physical particles in the bunch
      - rms_bunch_size: RMS bunch size at t=0 (vector) [m]
      - rms_velocity=[0,0,0]: RMS velocity spread at t=0 (vector) [m/s]
      - centroid_position=[0,0,0]: Position of the bunch centroid at t=0 (vector) [m]
      - centroid_velocity=[0,0,0]: Velocity (gamma*V) of the bunch centroid at t=0 (vector) [m/s]
      - velocity_divergence=[0,0,0]: Expansion rate of the bunch at t=0 (vector) [m/s/m]
    """
    def __init__(self,n_physical_particles, rms_bunch_size,
                 rms_velocity = [0.,0.,0.],
                 centroid_position = [0.,0.,0.],
                 centroid_velocity = [0.,0.,0.],
                 velocity_divergence = [0.,0.,0.],
                 **kw):
        self.n_physical_particles = n_physical_particles
        self.rms_bunch_size = rms_bunch_size
        self.rms_velocity = rms_velocity
        self.centroid_position = centroid_position
        self.centroid_velocity = centroid_velocity
        self.velocity_divergence = velocity_divergence

        self.init(**kw)


class PICMI_UniformDistribution(_ClassWithInit):
    """
    Describes a uniform density distribution of particles
      - density: Physical number density [m^-3]
      - lower_bound=[None,None,None]: Lower bound of the distribution (vector) [m]
      - upper_bound=[None,None,None]: Upper bound of the distribution (vector) [m]
      - rms_velocity=[0,0,0]: Thermal velocity spread (vector) [m/s]
      - directed_velocity=[0,0,0]: Directed, average, velocity (vector) [m/s]
      - fill_in=False: Flags whether to fill in the empty spaced opened up when the grid moves
    """

    def __init__(self, density,
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 fill_in = False,
                 **kw):
        self.density = density
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.rms_velocity = rms_velocity
        self.directed_velocity = directed_velocity
        self.fill_in = fill_in

        self.init(**kw)


class PICMI_AnalyticDistribution(_ClassWithInit):
    """
    Describes a uniform density plasma
      - density_expression: Analytic expression describing physical number density (string) [m^-3]
                            Expression should be in terms of the position, written as 'x', 'y', and 'z'.
      - lower_bound=[None,None,None]: Lower bound of the distribution (vector) [m]
      - upper_bound=[None,None,None]: Upper bound of the distribution (vector) [m]
      - rms_velocity=[0,0,0]: Thermal velocity spread (vector) [m/s]
      - directed_velocity=[0,0,0]: Directed, average, velocity (vector) [m/s]
      - fill_in=False: Flags whether to fill in the empty spaced opened up when the grid moves
    """

    def __init__(self, density_expression,
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 fill_in = False,
                 **kw):
        self.density_expression = density_expression
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.rms_velocity = rms_velocity
        self.directed_velocity = directed_velocity
        self.fill_in = fill_in

        self.init(**kw)


class PICMI_ParticleList(_ClassWithInit):
    """
    Load particles at the specified positions and velocities
      - x=0.: List of x positions of the particles [m]
      - y=0.: List of y positions of the particles [m]
      - z=0.: List of z positions of the particles [m]
      - ux=0.: List of ux positions of the particles (ux = gamma*vx) [m/s]
      - uy=0.: List of uy positions of the particles (uy = gamma*vy) [m/s]
      - uz=0.: List of uz positions of the particles (uz = gamma*vz) [m/s]
      - weight: Particle weight or list of weights, number of real particles per simulation particle
    """
    def __init__(self, x=0., y=0., z=0., ux=0., uy=0., uz=0., weight=0.,
                 **kw):
        # --- Get length of arrays, set to one for scalars
        lenx = np.size(x)
        leny = np.size(y)
        lenz = np.size(z)
        lenux = np.size(ux)
        lenuy = np.size(uy)
        lenuz = np.size(uz)
        lenw = np.size(weight)

        maxlen = max(lenx, leny, lenz, lenux, lenuy, lenuz, lenw)
        assert lenx==maxlen or lenx==1, "Length of x doesn't match len of others"
        assert leny==maxlen or leny==1, "Length of y doesn't match len of others"
        assert lenz==maxlen or lenz==1, "Length of z doesn't match len of others"
        assert lenux==maxlen or lenux==1, "Length of ux doesn't match len of others"
        assert lenuy==maxlen or lenuy==1, "Length of uy doesn't match len of others"
        assert lenuz==maxlen or lenuz==1, "Length of uz doesn't match len of others"
        assert lenw==maxlen or lenw==1, "Length of weight doesn't match len of others"

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
        # --- Note that weight can be a scalar

        self.weight = weight
        self.x = x
        self.y = y
        self.z = z
        self.ux = ux
        self.uy = uy
        self.uz = uz

        self.init(**kw)


# ------------------
# Numeric Objects
# ------------------


class PICMI_ParticleDistributionPlanarInjector(_ClassWithInit):
    """
    Describes the injection of particles from a plane
      - position: Position of the particle centroid (vector) [m]
      - plane_normal: Vector normal to the plane of injection (vector) [1]
      - plane_velocity: Velocity of the plane of injection (vector) [m/s]
      - method: InPlace - method of injection. One of 'InPlace', or 'Plane'
    """
    def __init__(self, position, plane_normal, plane_velocity=[0.,0.,0.], method='InPlace', **kw):
        self.position = position
        self.plane_normal = plane_normal
        self.plane_velocity = plane_velocity
        self.method = method

        self.init(**kw)


class PICMI_GriddedLayout(_ClassWithInit):
    """
    Specifies a gridded layout of particles
    - grid: grid object specifying the grid to follow
    - n_macroparticle_per_cell: number of particles per cell along each axis (vector)
    """
    def __init__(self, grid, n_macroparticle_per_cell, **kw):
        self.grid = grid
        self.n_macroparticle_per_cell = n_macroparticle_per_cell

        self.init(**kw)


class PICMI_PseudoRandomLayout(_ClassWithInit):
    """
    Specifies a pseudo-random layout of the particles

    Only one of these should be specified:
    - n_macroparticles: total number of macroparticles to load
    - n_macroparticles_per_cell: number of macroparticles to load per cell

    Optional argument
    - seed: pseudo-random number generator seed
    """
    def __init__(self, n_macroparticles=None, n_macroparticles_per_cell=None, seed=None, **kw):

        assert (n_macroparticles is not None)^(n_macroparticles_per_cell is not None), \
               Exception('Only one of n_macroparticles and n_macroparticles_per_cell must be specified')

        self.n_macroparticles = n_macroparticles
        self.n_macroparticles_per_cell = n_macroparticles_per_cell
        self.seed = seed

        self.init(**kw)


class PICMI_BinomialSmoother(_ClassWithInit):
    """
    Descibes a binomial smoother operator (applied to grids)
    - n_pass: Number of passes along each axis (vector)
    - compensation: Flags whether to apply comensation (logical)
    """
    def __init__(self, n_pass=None, compensation=None, **kw):
        self.n_pass = n_pass
        self.compensation = compensation

        self.init(**kw)


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

        self.init(**kw)


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

        self.init(**kw)


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

        self.init(**kw)


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

    Methods:
      - add: Add object to solver (e.g. laser)
    """

    methods_list = ['Yee', 'CKC', 'Lehe', 'PSTD', 'PSATD', 'GPSTD']

    def __init__(self, grid, method=None, stencil_order=None, cfl=None, l_nodal=None,
                 source_smoother=None, field_smoother=None,
                 **kw):

        assert method is None or method in PICMI_ElectromagneticSolver.Methods_list, \
               Exception('method must be one of '+', '.join(PICMI_ElectromagneticSolver.methods_list))

        self.grid = grid
        self.method = method
        self.cfl = cfl
        self.stencil_order = stencil_order
        self.l_nodal = l_nodal
        self.source_smoother = source_smoother
        self.field_smoother = field_smoother

        self.init(**kw)


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

        self.init(**kw)


# Laser related objects


class PICMI_GaussianLaser(_ClassWithInit):
    """
    Specifies a Gaussian laser distribution
      - wavelength: Laser wavelength
      - waist: Waist of the Gaussian pulse at focus [m]
      - duration: Duration of the Gaussian pulse [s]
      - focal_position=[0,0,0]: Position of the laser focus (vector) [m]
      - centroid_position=[0,0,0]: Position of the laser centroid at time 0 (vector) [m]
      - propagation_direction=[0,0,1]: Direction of propagation (unit vector) [1]
      - polarization_angle=0: Angle of polarization (relative to X) [radians]
      - a0: Normalized vector potential at focus
            Specify either a0 or E0 (E0 takes precedence).
      - E0: Maximum amplitude of the laser field [V/m]
            Specify either a0 or E0 (E0 takes precedence).
    """
    def __init__(self, wavelength, waist, duration,
                 focal_position = [0., 0., 0.],
                 centroid_position = [0., 0., 0.],
                 propagation_direction = [0., 0., 1.],
                 polarization_angle = 0.,
                 a0 = None, E0 = None,
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
        self.centroid_position = centroid_position
        self.propagation_direction = propagation_direction
        self.polarization_angle = polarization_angle
        self.a0 = a0
        self.E0 = E0

        self.init(**kw)


class PICMI_LaserAntenna(_ClassWithInit):
    """
    Specifies the laser antenna injection method
      - position: Position of antenna launching the laser (vector) [m]
      - normal_vector: Vector normal to antenna plane (vector) [1]
    """
    def __init__(self, position, normal_vector, **kw):

        self.position = position
        self.normal_vector = normal_vector

        self.init(**kw)


# Main simuation object


class PICMI_Simulation(_ClassWithInit):
    """
    Simulation
      - solver: Field solver to be used in the simulation (an instance of one of the implemented solvers)
      - time_step_size: Absolute time step size of the simulation [s]
                        (needed if the CFL is not specified elsewhere)
      - max_steps: Maximum number of time steps
      - max_time: Maximum time to run the simulation [s]
      - verbose: Verbosity flag
      - particle_shape: Default particle shape for species added to this simulation. Possible values are 'NGP', 'linear', 'quadratic', 'cubic'
    """

    def __init__(self, solver=None, time_step_size=None, max_steps=None, max_time=None, verbose=None,
                particle_shape='linear', **kw):

        self.solver = solver
        self.time_step_size = time_step_size
        self.verbose = verbose
        self.max_steps = max_steps
        self.max_time = max_time

        self.species = []
        self.layouts = []
        self.initialize_self_fields = []

        self.lasers = []
        self.laser_injection_methods = []

        self.init(**kw)

    def add_species(self, species, layout, initialize_self_field=False):
        """
        Add species to be used in the simulation
        - species: species object
        - layout: particle layout for initial distribution
        - initialize_self_field=True: whether the initial space-charge fields
        of this species is calculated and added to the simulation.
        """
        self.species.append(species)
        self.layouts.append(layout)
        self.initialize_self_fields.append(initialize_self_field)

    def add_laser(self, laser, injection_method):
        """
        Add a laser pulses that to be injected in the simulation
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
        self.lasers.append(laser)
        self.laser_injection_methods.append(injection_method)

    def write_input_file(self, file_name):
        raise NotImplementedError

    def step(self, nsteps=1):
        raise NotImplementedError
