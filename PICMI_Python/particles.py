"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
The classes in the file are all particle related
"""
import math
import sys
import re
import numpy as np

from .base import _ClassWithInit

# ---------------
# Physics objects
# ---------------


class PICMI_Species(_ClassWithInit):
    """
    Sets up the species to be simulated.
    The species charge and mass can be specified by setting the particle type or by setting them directly.
    If the particle type is specified, the charge or mass can be set to override the value from the type.

    Parameters
    ----------
    particle_type: string, optional
        A string specifying an elementary particle, atom, or other, as defined in
        the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md

    name: string, optional
        Name of the species

    method: {'Boris', 'Vay', 'Higuera-Cary', 'Li' , 'free-streaming', 'LLRK4'}
        The particle advance method to use. Code-specific method can be specified using 'other:<method>'. The default is code
        dependent.

        - 'Boris': Standard "leap-frog" Boris advance
        - 'Vay':
        - 'Higuera-Cary':
        - 'Li' :
        - 'free-streaming': Advance with no fields
        - 'LLRK4': Landau-Lifschitz radiation reaction formula with RK-4)

    charge_state: float, optional
        Charge state of the species (applies only to atoms) [1]

    charge: float, optional
        Particle charge, required when type is not specified, otherwise determined from type [C]

    mass: float, optional
        Particle mass, required when type is not specified, otherwise determined from type [kg]

    initial_distribution: distribution instance
        The initial distribution loaded at t=0. Must be one of the standard distributions implemented.

    density_scale: float, optional
        A scale factor on the density given by the initial_distribution

    particle_shape: {'NGP', 'linear', 'quadratic', 'cubic'}
        Particle shape used for deposition and gather.
        If not specified, the value from the `Simulation` object will be used.
        Other values maybe specified that are code dependent.
    """

    methods_list = ['Boris' , 'Vay', 'Higuera-Cary', 'Li', 'free-streaming', 'LLRK4']

    def __init__(self, particle_type=None, name=None, charge_state=None, charge=None, mass=None,
                 initial_distribution=None, particle_shape=None, density_scale=None, method=None, **kw):


        assert method is None or method in PICMI_Species.methods_list or method.startswith('other:'), \
            Exception('method must starts with either "other:", or be one of the following '+', '.join(PICMI_Species.methods_list))

        self.method = method
        self.particle_type = particle_type
        self.name = name
        self.charge = charge
        self.charge_state = charge_state
        self.mass = mass
        self.initial_distribution = initial_distribution
        self.particle_shape = particle_shape
        self.density_scale = density_scale

        self.interactions = []

        self.handle_init(kw)


class PICMI_MultiSpecies(_ClassWithInit):
    """
    INCOMPLETE: proportions argument is not implemented
    Multiple species that are initialized with the same distribution.
    Each parameter can be list, giving a value for each species, or a single value which is given to all species.
    The species charge and mass can be specified by setting the particle type or by setting them directly.
    If the particle type is specified, the charge or mass can be set to override the value from the type.

    Parameters
    ----------
    particle_types: list of strings, optional
        A string specifying an elementary particle, atom, or other, as defined in
        the openPMD 2 species type extension, openPMD-standard/EXT_SpeciesType.md

    names: list of strings, optional
        Names of the species

    charge_states: list of floats, optional
        Charge states of the species (applies only to atoms)

    charges: list of floats, optional
        Particle charges, required when type is not specified, otherwise determined from type [C]

    masses: list of floats, optional
        Particle masses, required when type is not specified, otherwise determined from type [kg]

    proportions: list of floats, optional
        Proportions of the initial distribution made up by each species

    initial_distribution: distribution instance
        Initial particle distribution, applied to all species

    particle_shape: {'NGP', 'linear', 'quadratic', 'cubic'}
        Particle shape used for deposition and gather.
        If not specified, the value from the `Simulation` object will be used.
        Other values maybe specified that are code dependent.
    """
    # --- Note to developer: This class attribute needs to be set to the Species class
    # --- defined in the codes PICMI implementation.
    Species_class = None

    def __init__(self, particle_types=None, names=None, charge_states=None, charges=None, masses=None,
                 proportions=None, initial_distribution=None, particle_shape=None,
                 **kw):

        self.particle_types = particle_types
        self.names = names
        self.charges = charges
        self.charge_states = charge_states
        self.masses = masses
        self.proportions = proportions
        self.initial_distribution = initial_distribution
        self.particle_shape = particle_shape

        self.nspecies = None
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
            particle_type = self.get_input_item(particle_types, i)
            name = self.get_input_item(names, i)
            charge = self.get_input_item(charges, i)
            charge_state = self.get_input_item(charge_states, i)
            mass = self.get_input_item(masses, i)
            proportion = self.get_input_item(proportions, i)
            specie = PICMI_MultiSpecies.Species_class(particle_type = particle_type,
                                                      name = name,
                                                      charge = charge,
                                                      charge_state = charge_state,
                                                      mass = mass,
                                                      initial_distribution = initial_distribution,
                                                      density_scale = proportion)
            self.species_instances_list.append(specie)
            if name is not None:
                self.species_instances_dict[name] = specie

        self.handle_init(kw)

    def check_nspecies(self, var):
        if var is not None:
            try:
                nvars = len(var)
            except TypeError:
                nvars = 1
            assert self.nspecies is None or self.nspecies == nvars, Exception('All inputs must have the same length')
            self.nspecies = nvars

    def get_input_item(self, var, i):
        if var is None:
            return None
        else:
            try:
                len(var)
            except TypeError:
                return var
            else:
                return var[i]

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

    Parameters
    ----------
    n_physical_particles: integer
        Number of physical particles in the bunch

    rms_bunch_size: vector of floats
        RMS bunch size at t=0 [m]

    rms_velocity: vector of floats, default=[0.,0.,0.]
        RMS velocity spread at t=0 [m/s]

    centroid_position: vector of floats, default=[0.,0.,0.]
        Position of the bunch centroid at t=0 [m]

    centroid_velocity: vector of floats, default=[0.,0.,0.]
        Velocity (gamma*V) of the bunch centroid at t=0 [m/s]

    velocity_divergence: vector of floats, default=[0.,0.,0.]
        Expansion rate of the bunch at t=0 [m/s/m]
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

        self.handle_init(kw)


class PICMI_UniformDistribution(_ClassWithInit):
    """
    Describes a uniform density distribution of particles

    Parameters
    ----------
    density: float
        Physical number density [m^-3]

    lower_bound: vector of floats, optional
        Lower bound of the distribution [m]

    upper_bound: vector of floats, optional
        Upper bound of the distribution [m]

    rms_velocity: vector of floats, default=[0.,0.,0.]
        Thermal velocity spread [m/s]

    directed_velocity: vector of floats, default=[0.,0.,0.]
        Directed, average, velocity [m/s]

    fill_in: bool, optional
        Flags whether to fill in the empty spaced opened up when the grid moves
    """

    def __init__(self, density,
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 fill_in = None,
                 **kw):
        self.density = density
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.rms_velocity = rms_velocity
        self.directed_velocity = directed_velocity
        self.fill_in = fill_in

        self.handle_init(kw)


class PICMI_UniformFluxDistribution(_ClassWithInit):
    """
    Describes a flux of particles emitted from a plane

    Parameters
    ----------
    flux: float
        Flux of particles [m^-2.s^-1]

    flux_normal_axis: string
        x, y, or z for 3D, x or z for 2D, or r, t, or z in RZ geometry

    surface_flux_position: double
        location of the injection plane [m] along the direction
        specified by `flux_normal_axis`

    flux_direction: int
        Direction of the flux relative to the plane: -1 or +1

    lower_bound: vector of floats, optional
        Lower bound of the distribution [m]

    upper_bound: vector of floats, optional
        Upper bound of the distribution [m]

    rms_velocity: vector of floats, default=[0.,0.,0.]
        Thermal velocity spread [m/s]

    directed_velocity: vector of floats, default=[0.,0.,0.]
        Directed, average, velocity [m/s]
    """

    def __init__(self, flux, flux_normal_axis,
                 surface_flux_position, flux_direction,
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 **kw):
        self.flux = flux
        self.flux_normal_axis = flux_normal_axis
        self.surface_flux_position = surface_flux_position
        self.flux_direction = flux_direction
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.rms_velocity = rms_velocity
        self.directed_velocity = directed_velocity

        self.handle_init(kw)

class PICMI_AnalyticDistribution(_ClassWithInit):
    """
    Describes a uniform density plasma

    Parameters
    ----------
    density_expression: string
        Analytic expression describing physical number density (string) [m^-3].
        Expression should be in terms of the position, written as 'x', 'y', and 'z'.
        Parameters can be used in the expression with the values given as keyword arguments.

    momentum_expressions: list of strings
        Analytic expressions describing the gamma*velocity for each axis [m/s].
        Expressions should be in terms of the position, written as 'x', 'y', and 'z'.
        Parameters can be used in the expression with the values given as keyword arguments.
        For any axis not supplied (set to None), directed_velocity will be used.

    lower_bound: vector of floats, optional
        Lower bound of the distribution [m]

    upper_bound: vector of floats, optional
        Upper bound of the distribution [m]

    rms_velocity: vector of floats, detault=[0.,0.,0.]
        Thermal velocity spread [m/s]

    directed_velocity: vector of floats, detault=[0.,0.,0.]
        Directed, average, velocity [m/s]

    fill_in: bool, optional
        Flags whether to fill in the empty spaced opened up when the grid moves


    This will create a distribution where the density is n0 below rmax and zero elsewhere.::

    .. code-block: python

      dist = AnalyticDistribution(density_expression='((x**2+y**2)<rmax**2)*n0',
                                  rmax = 1.,
                                  n0 = 1.e20,
                                  ...)

    """

    def __init__(self, density_expression,
                 momentum_expressions = [None, None, None],
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 fill_in = None,
                 **kw):
        self.density_expression = '{}'.format(density_expression).replace('\n', '')
        self.momentum_expressions = momentum_expressions
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.rms_velocity = rms_velocity
        self.directed_velocity = directed_velocity
        self.fill_in = fill_in

        # --- Convert momentum_expressions to string if needed.
        if self.momentum_expressions[0] is not None:
            self.momentum_expressions[0] = '{}'.format(self.momentum_expressions[0]).replace('\n', '')
        if self.momentum_expressions[1] is not None:
            self.momentum_expressions[1] = '{}'.format(self.momentum_expressions[1]).replace('\n', '')
        if self.momentum_expressions[2] is not None:
            self.momentum_expressions[2] = '{}'.format(self.momentum_expressions[2]).replace('\n', '')

        # --- Find any user defined keywords in the kw dictionary.
        # --- Save them and delete them from kw.
        # --- It's up to the code to make sure that all parameters
        # --- used in the expression are defined.
        self.user_defined_kw = {}
        for k in list(kw.keys()):
            if re.search(r'\b%s\b'%k, self.density_expression):
                self.user_defined_kw[k] = kw[k]
                del kw[k]
            elif self.momentum_expressions[0] is not None and re.search(r'\b%s\b'%k, self.momentum_expressions[0]):
                self.user_defined_kw[k] = kw[k]
                del kw[k]
            elif self.momentum_expressions[1] is not None and re.search(r'\b%s\b'%k, self.momentum_expressions[1]):
                self.user_defined_kw[k] = kw[k]
                del kw[k]
            elif self.momentum_expressions[2] is not None and re.search(r'\b%s\b'%k, self.momentum_expressions[2]):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

        self.handle_init(kw)


class PICMI_ParticleListDistribution(_ClassWithInit):
    """
    Load particles at the specified positions and velocities

    Parameters
    ----------
    x: float, default=0.
        List of x positions of the particles [m]

    y: float, default=0.
        List of y positions of the particles [m]

    z: float, default=0.
        List of z positions of the particles [m]

    ux: float, default=0.
        List of ux positions of the particles (ux = gamma*vx) [m/s]

    uy: float, default=0.
        List of uy positions of the particles (uy = gamma*vy) [m/s]

    uz: float, default=0.
        List of uz positions of the particles (uz = gamma*vz) [m/s]

    weight: float
        Particle weight or list of weights, number of real particles per simulation particle
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

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


class PICMI_ParticleDistributionPlanarInjector(_ClassWithInit):
    """
    Describes the injection of particles from a plane

    Parameters
    ----------
    position: vector of floats
        Position of the particle centroid [m]

    plane_normal: vector of floats
        Vector normal to the plane of injection [1]

    plane_velocity: vector of floats
        Velocity of the plane of injection [m/s]

    method: {'InPlace', 'Plane'}
    """
    def __init__(self, position, plane_normal, plane_velocity=[0.,0.,0.], method='InPlace', **kw):
        self.position = position
        self.plane_normal = plane_normal
        self.plane_velocity = plane_velocity
        self.method = method

        self.handle_init(kw)


class PICMI_GriddedLayout(_ClassWithInit):
    """
    Specifies a gridded layout of particles

    Parameters
    ----------
    n_macroparticle_per_cell: vector of integers
        Number of particles per cell along each axis

    grid: grid instance, optional
        Grid object specifying the grid to follow.
        If not specified, the underlying grid of the code is used.
    """
    def __init__(self, n_macroparticle_per_cell, grid=None, **kw):
        self.n_macroparticle_per_cell = n_macroparticle_per_cell
        self.grid = grid

        self.handle_init(kw)


class PICMI_PseudoRandomLayout(_ClassWithInit):
    """
    Specifies a pseudo-random layout of the particles

    Parameters
    ----------
    n_macroparticles: integer
        Total number of macroparticles to load.
        Either this argument or n_macroparticles_per_cell should be supplied.

    n_macroparticles_per_cell: integer
        Number of macroparticles to load per cell.
        Either this argument or n_macroparticles should be supplied.

    seed: integer, optional
        Pseudo-random number generator seed

    grid: grid instance, optional
        Grid object specifying the grid to follow for n_macroparticles_per_cell.
        If not specified, the underlying grid of the code is used.
    """
    def __init__(self, n_macroparticles=None, n_macroparticles_per_cell=None, seed=None, grid=None, **kw):

        assert (n_macroparticles is not None)^(n_macroparticles_per_cell is not None), \
               Exception('Only one of n_macroparticles and n_macroparticles_per_cell must be specified')

        self.n_macroparticles = n_macroparticles
        self.n_macroparticles_per_cell = n_macroparticles_per_cell
        self.seed = seed
        self.grid = grid

        self.handle_init(kw)
