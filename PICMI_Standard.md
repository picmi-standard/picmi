The Particle-In-Cell Modeling Interface (PICMI) Standard
==================================================

VERSION: **0.0.2** (November 22th, 2017)

Code specific specifications
----------------------------
  - `codename`
    - **type**: *string*
    - **definition**: "Name of code"

Physics objects
---------------

### Particle species and distributions

  - `Species`
    - **type**: *object*
    - `type` - **type**: *string* - "an elementary particle or atom, or user-defined type, as defined in the openPMD species type extension"
    - `name` - **type**: *string* - "name of the species."
    - `charg_state` (optional) - **type**: *integer* - "Charge state of the species (applies to atoms)."
    - `charge` (optional) - **type**: *double* - **default**: Type.Charge (multiplied by Charge_state if atom) - "particles charge, overwrites default if provided."
    - `mass` (optional) - **type**: *double* - **default**: Type.Mass - "particles mass, overwrites default if provided."

  - `GaussianBeam`
    - **type**: *object*
    - `species` - **type**: *Particle* - "Particle species"
    - `number_real_particles` - **type**: *double* - "Number of real particles in the beam."
    - `number_sim_particles` - **type**: *double* - "Number of simulation particles in the beam."
    - `T0` - **type**: *double* - **default**: 0. - "Time at which parameters are specified [s]."
    - `Xmean` - **type**: *double* - **default**: 0. - "Mean X position [m]."
    - `Ymean` - **type**: *double* - **default**: 0. - "Mean Y position [m]."
    - `Zmean` - **type**: *double* - **default**: 0. - "Mean Z position [m]."
    - `Xrms` - **type**: *double* - **default**: 0. - "R.M.S. size along X [m]."
    - `Yrms` - **type**: *double* - **default**: 0. - "R.M.S. size along Y [m]."
    - `Zrms` - **type**: *double* - **default**: 0. - "R.M.S. size along Z [m]."
    - `UXmean` - **type**: *double* - **default**: 0. - "Mean velocity (gamma*V) along X [m/s]."
    - `UYmean` - **type**: *double* - **default**: 0. - "Mean velocity (gamma*V) along X [m/s]."
    - `UZmean` - **type**: *double* - **default**: 0. - "Mean velocity (gamma*V) along X [m/s]."
    - `UXrms` - **type**: *double* - **default**: 0. - "R.M.S. velocity (gamma*V) spread along X [m/s]."
    - `UYrms` - **type**: *double* - **default**: 0. - "R.M.S. velocity (gamma*V) spread along Y [m/s]."
    - `UZrms` - **type**: *double* - **default**: 0. - "R.M.S. velocity (gamma*V) spread along Z [m/s]."
    - `UXdiv` - **type**: *double* - **default**: 0. - "Velocity (gamma*V) divergence along X [m/s/m]."
    - `UYdiv` - **type**: *double* - **default**: 0. - "Velocity (gamma*V) divergence along Y [m/s/m]."
    - `UZdiv` - **type**: *double* - **default**: 0. - "Velocity (gamma*V) divergence along Z [m/s/m]."
    - `density_func` (optional) - **type**: *function*: "Function modulating density as a function of x, y, z and/or time."
    - `array_func` (optional) - **type**: *array*: "Array modulating density as a function of x, y, z and/or time."

  - `Plasma`
    - **type**: *object*
    - `species` - **type**: *Particle list* - "Particle species or list of species"
    - `density` - **type**: *double* - "Plasma density [m^-3]."
    - `xmin` - **type**: *double* - **default**: -infinity - "Min position of box along X."
    - `xmax` - **type**: *double* - **default**: +infinity - "Max position of box along X."
    - `ymin` - **type**: *double* - **default**: -infinity - "Min position of box along Y."
    - `ymax` - **type**: *double* - **default**: +infinity - "Max position of box along Y."
    - `zmin` - **type**: *double* - **default**: -infinity - "Min position of box along Z."
    - `zmax` - **type**: *double* - **default**: +infinity - "Max position of box along Z."
    - `vthx` - **type**: *double* - **default**: 0. - "Thermal velocity along X."
    - `vthy` - **type**: *double* - **default**: 0. - "Thermal velocity along Y."
    - `vthz` - **type**: *double* - **default**: 0. - "Thermal velocity along Z."
    - `vxmean` - **type**: *double* - **default**: 0. - "Mean velocity along X."
    - `vymean` - **type**: *double* - **default**: 0. - "Mean velocity along Y."
    - `vzmean` - **type**: *double* - **default**: 0. - "Mean velocity along Z."
    - `number_per_cell` - **type**: *double* - "Number of particles per cell (randomly placed)
                                                Only one of number_per_cell or number_per_cell_each_dim should be specified."
    - `number_per_cell_each_dim` - **type**: *integer array* - **size**: Ndims - "Number of particles along each axis (for regularly placed particles)
                                                Only one of number_per_cell or number_per_cell_each_dim should be specified."
    - `density_func` (optional) - **type**: *function*: "Function modulating density as a function of x, y, z and/or time."
    - `array_func` (optional) - **type**: *array*: "Array modulating density as a function of x, y, z and/or time."

  - `ParticleList`
    - **type**: *object*
    - `species` - **type**: *Particle* - "Particle species"
    - `weight` - **type**: *double* - "Particle weight, number of real particles per simulation particle"
    - `x` - **type**: *list* - "List of x positions of the particles [m]"
    - `y` - **type**: *list* - "List of y positions of the particles [m]"
    - `z` - **type**: *list* - "List of z positions of the particles [m]"
    - `ux` - **type**: *list* - "List of ux positions of the particles (ux = gamma*vx) [m/s]"
    - `uy` - **type**: *list* - "List of uy positions of the particles (uy = gamma*vy) [m/s]"
    - `uz` - **type**: *list* - "List of uz positions of the particles (uz = gamma*vz) [m/s]"


Numerics objects
---------------
### Particles

  - `ParticleDistributionInjector`
    - **type**: *object*
    - `distribution` - **type**: *Beam or Plasma* - "beam/plasma to inject"
    - `method` - **type**: *string* - **default**: "InPlace" - "method of injection ('InPlace', 'Plane')"
    - `X0` - **type**: *double* - **default**: 0. - "Position of the particle centroid in X."
    - `Y0` - **type**: *double* - **default**: 0. - "Position of the particle centroid in Y."
    - `Z0` - **type**: *double* - **default**: 0. - "Position of the particle centroid in Z."
    - `Xplane` (optional) - **type**: *double* - **default**: 0. - "Position of the plane of injection in X."
    - `Yplane` (optional) - **type**: *double* - **default**: 0. - "Position of the plane of injection in Y."
    - `Zplane` (optional) - **type**: *double* - **default**: 0. - "Position of the plane of injection in Z."
    - `VXplane` (optional) - **type**: *double* - **default**: 0. - "Velocity of the plane of injection in X."
    - `VYplane` (optional) - **type**: *double* - **default**: 0. - "Velocity of the plane of injection in Y."
    - `VZplane` (optional) - **type**: *double* - **default**: 0. - "Velocity of the plane of injection in Z."
    - `XVecPlane` - **type**: *double* - "Component along X of vector normal to injection plane."
    - `YVecPlane` - **type**: *double* - "Component along Y of vector normal to injection plane."
    - `ZVecPlane` - **type**: *double* - "Component along Z of vector normal to injection plane."

###Fields
  - `Grid`
    - **type**: *object*
    - `nx` - **type**: *integer* - "Number of cells along X (Nb nodes=nx+1)."
    - `ny` - **type**: *integer* - "Number of cells along Y (Nb nodes=ny+1)."
    - `nr` - **type**: *integer* - "Number of cells along R (Nb nodes=nr+1)."
    - `nz` - **type**: *integer* - "Number of cells along Z (Nb nodes=nz+1)."
    - `nm` - **type**: *integer* - "Number of azimuthal modes."
    - `xmin` - **type**: *double* - "Position of first node along X."
    - `xmax` - **type**: *double* - "Position of last node along X."
    - `ymin` - **type**: *double* - "Position of first node along Y."
    - `ymax` - **type**: *double* - "Position of last node along Y."
    - `rmax` - **type**: *double* - "Position of last node along R."
    - `zmin` - **type**: *double* - "Position of first node along Z."
    - `zmax` - **type**: *double* - "Position of last node along Z."
    - `bcxmin` - **type**: *string* - "Boundary condition at min X: periodic/open/dirichlet/neumann."
    - `bcxmax` - **type**: *string* - "Boundary condition at max X: periodic/open/dirichlet/neumann."
    - `bcymin` - **type**: *string* - "Boundary condition at min Y: periodic/open/dirichlet/neumann."
    - `bcymax` - **type**: *string* - "Boundary condition at max Y: periodic/open/dirichlet/neumann."
    - `bcrmax` - **type**: *string* - "Boundary condition at max R: open/dirichlet/neumann."
    - `bczmin` - **type**: *string* - "Boundary condition at min Z: periodic/open/dirichlet/neumann."
    - `bczmax` - **type**: *string* - "Boundary condition at max Z: periodic/open/dirichlet/neumann."
    - `moving_window_velocity` - **type** *double array* - **size**: Ndims - "An array of the moving frame velocity in each direction"

    - `getmins()` - **type**: *method* - "Get mins of grid"
    - `getmaxs()` - **type**: *method* - "Get maxs of grid"

  - `EM_solver`
    - **type**: *Grid*
    - `method` - **type**: *string* - "Yee/CK/CKC/Lehe/PSTD/PSATD/GPSTD"
    - `norderx` - **type**: *integer* - "Order of stencil in X (-1=infinite)."
    - `nordery` - **type**: *integer* - "Order of stencil in Y (-1=infinite)."
    - `norderr` - **type**: *integer* - "Order of stencil in R (-1=infinite)."
    - `norderz` - **type**: *integer* - "Order of stencil in Z (-1=infinite)."
    - `l_nodal` - **type**: *logical* - "Quantities are at nodes if True, staggered otherwise."
    - `add` - **type**: *method* - "Add object to solver (e.g. laser)"

  - `ES_solver`
    - **type**: *Grid*
    - `method` - **type**: *double* - "FFT/Multigrid"

### Simulation
  - `Simulation`
    - **Type**: *object*
    - **Input arguments:**
        - `timestep` - **type** *float* - **default**: 0. - "Absolute time step size of the simulation
        (use 0 if you prefer specifying instead the timestep relative to the CFL limit)"
        - `timestep_over_CFL` - **type** *float* - **default**: 1. - "Ratio of the time step size to the CFL limit
        (used only if `timestep` is 0 ; should raise an error when the code does not have a well-defined CFL)"
        - `max_step` - **type** *integer* - "Maximum number of time steps"
        - `max_time` - **type** *float* - "Maximum time to run the simulation"
        - `verbose` - **type** *boolean* - "Verbosity flag"
    - **Methods:**
        - `step(`
        `nsteps` - **type** *integer* - "Number of time steps"
        `)`

Injecting a laser pulse in the Simulation
-----------------------------------------

Laser pulses are injected in the simulation by using the function `add_laser_pulse`:

  - `add_laser_pulse( simulation, laser_profile, injection_method )`:
    - `simulation`: a `Simulation` object (see below)
        Top-level object that contains all the relevant data for a simulation.
    - `laser_profile`: one of laser profile object (see below)
        Specifies the **physical** properties of the laser pulse.
        (e.g. spatial and temporal profile, wavelength, amplitude, etc.)
    - `injection_method`: a laser injector object, optional (see below)
        Specifies how the laser is injected (numerically) into the simulation
        (e.g. through a laser antenna, or directly added to the mesh).
        This argument describes an **algorithm**, not a physical object.
        It is optional. (It is up to each code to define the default method
        of injection, if the user does not provide `injection_method`)

### Defining the physical parameters, through a laser profile object

Instances of the different classes below can be passed as the `laser_profile`
argument in `add_laser_pulse`:

  - `GaussianLaser`
    - `wavelength` - **type**: *double* - "Laser wavelength."
    - `waist` - **type**: *double* - "Waist of the Gaussian pulse at focus [m]."
    - `duration` - **type**: *double* - "Duration of the Gaussian pulse [s]."
    - `focal_position` - **type**: *double* - "Position of the laser focus."
    - `x0` - **type**: *double* - "Position of the laser centroid in X at time 0."
    - `y0` - **type**: *double* - "Position of the laser centroid in Y at time 0."
    - `z0` - **type**: *double* - "Position of the laser centroid in Z at time 0"
    - `pol_angle` - **type**: *double* - "Angle of polarization (relative to X)."
    - `a0` - **type**: *double* - "Normalized vector potential at focus. Specify eiter a0 or E0."
    - `E0` - **type**: *double* - "Maximum amplitude of the laser field (in V/m). Specify eiter a0 or E0."

### Defining an injection method, through a laser injector object

Instances of the different classes below can be passed as the `injection_method`
argument of `add_laser_pulse`:

  - `LaserAntenna`
    - `laser` - **type**: *Laser* - "Laser to be injected."
    - `antenna_x0` - **type**: *double* - "Position of antenna launching the laser along X."
    - `antenna_y0` - **type**: *double* - "Position of antenna launching the laser along Y."
    - `antenna_z0` - **type**: *double* - "Position of antenna launching the laser along Z."
    - `antenna_xvec` - **type**: *double* - "Component along X of vector normal to antenna plane."
    - `antenna_yvec` - **type**: *double* - "Component along Y of vector normal to antenna plane."
    - `antenna_zvec` - **type**: *double* - "Component along Z of vector normal to antenna plane."

Examples
--------

### Example 1
Beam_Electrons = Species(type='electron')

EBeam = GaussianBeam(species=Beam_Electrons,
                     nbpart=1.e9,
                     Xrms=1.e-6,
                     Yrms=1.e-6,
                     Zrms=5.e-6,
                     UXrms=1.e3,
                     UYrms=1.e3,
                     UZrms=1.e3,
                     UZmean=1.e8)

LBeam = GaussianLaser(wavelength=1.e-6,
                      waist=30.e-6,
                      duration=30.e-15,
                      pol_angle=pi/2,
                      a0=2.)

EBeam_Injector = ParticleDistributionInjector(distribution=EBeam,
                                             method='InPlace',
                                             Z0=10.e-6)

LBeam_Antenna = LaserAntenna(laser=LBeam,
                             focal_position=40.e-6,
                             z0=40.e-6-LBeam.duration*clight,
                             antenna_z0=30.e-6)

To define: grid, sim.
