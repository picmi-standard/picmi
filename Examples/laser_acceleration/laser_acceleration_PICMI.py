# This should be the only line that needs to be changed for different codes
# e.g. `from pywarpx import picmi`
#      `from fbpic import picmi`
#      `from warp import picmi`
from pywarpx import picmi

import picmi.constants as cst

# Run parameters - can be in separate file
# ========================================

# Physics parameters
# ------------------

# --- laser
laser_a0             = 4.        # Normalized potential vector
laser_wavelength     = 8e-07     # Wavelength of the laser (in meters)
laser_waist          = 5e-06     # Waist of the laser (in meters)
laser_duration       = 15e-15    # Duration of the laser (in seconds)
laser_polarization   = cst.pi/2.  # Polarization angle (in rad)
laser_injection_loc  = 9.e-6     # Position of injection (in meters, along z)
laser_focal_distance = 100.e-6   # Focal distance from the injection (in meters)
laser_t_peak         = 30.e-15   # The time at which the laser reaches its peak
                                  # at the antenna injection location (in seconds)
# --- plasma
plasma_density_expression = "1.e23*(1+tanh((z - 20.e-6)/10.e-6))/2."
plasma_min     = [-20.e-6, -20.e-6,  0.0e-6]
plasma_max     = [ 20.e-6,  20.e-6,  1.e-3]

# --- electron bunch
bunch_physical_particles  = 1.e9
bunch_rms_size            = [2.e-6, 2.e-6, 1.e-6]
bunch_rms_velocity        = [0.,0.,10.]
bunch_centroid_position   = [0.,0.,100.e-6]
bunch_centroid_velocity   = [0.,0.,1000.]
bunch_velocity_divergence = [0.,0.,0.]

# Numerics parameters
# -------------------

# --- Nb time steps
max_steps = 1000

# --- grid
nx = 64
ny = 64
nz = 480
xmin = 1.5*plasma_min[0]
xmax = 1.5*plasma_max[0]
ymin = 1.5*plasma_min[1]
ymax = 1.5*plasma_max[1]
zmin = -56.e-6
zmax = 12.e-6
moving_window_velocity = [0., 0., cst.c]
number_per_cell_each_dim = [2, 2, 1]

# --- geometry and solver
em_solver_method = 'CKC'
geometry = '3D'
# Note that code-specific change can be introduced with `picmi.codename`
if picmi.codename == 'fbpic':
    em_solver_method = 'PSATD'
    geometry = 'RZ'

# Physics part - can be in separate file
# ======================================

# Physics components
# ------------------

# --- laser

laser = picmi.GaussianLaser(
    wavelength            = laser_wavelength,
    waist                 = laser_waist,
    duration              = laser_duration,
    focal_position        = [0., 0., laser_focal_distance + laser_injection_loc],
    centroid_position     = [0., 0., laser_injection_loc - cst.c*laser_t_peak],
    polarization_angle    = laser_polarization,
    propagation_direction = [0,0,1],
    a0                    = laser_a0))

laser_antenna = picmi.LaserAntenna(position      = [0., 0., laser_injection_loc],  # This point is on the laser plane
                                   normal_vector = [0., 0., 1.])  # The plane normal direction

# --- plasma

uniform_plasma = picmi.UniformDistribution(density     = plasma_density,
                                           lower_bound = plasma_min,
                                           upper_bound = plasma_max,
                                           fill_in     = True)
plasma_dist = picmi.AnalyticDistribution(
                density_expression = plasma_density,
                fill_in=True)


plasma = picmi.MultiSpecies(
                particle_types = ['He', 'Ar', 'electron'],
                names          = ['He+', 'Argon', 'e-'],
                charge_states  = [1, 5, None],
                proportions    = [0.2, 0.8, 0.2 + 5*0.8],
                initial_distribution=plasma_dist)

# Individual species in a `MultiSpecies` can be addressed either
# with their index (using Python indexing conventions) or with their name
# (if the user provided a name)
# Set the ionization for the species number 1 (Argon)
# and place the created electrons into the species number 2 (electron)
if picmi.codename is not 'pywarpx':
    plasma['Argon'].activate_ionization(model="ADK", target_species=plasma['e-'])

# --- electron beam

beam = picmi.Species(
            particle_type        = 'electron',
            name                 = 'beam',
            initial_distribution = beam_dist)


# Numerics components
# -------------------

if geometry == '3D':
    grid = picmi.Cartesian3DGrid(number_of_cells = [nx, ny, nz],
                                 lower_bound = [xmin, ymin, zmin],
                                 upper_bound = [xmax, ymax, zmax],
                                 lower_boundary_conditions = ['periodic', 'periodic', 'open'],
                                 upper_boundary_conditions = ['periodic', 'periodic', 'open'],
                                 moving_window_velocity = moving_window_velocity,
                                 warpx_max_grid_size=32)
elif geometry == 'RZ':
    nr = nx/2.
    grid = picmi.CylindricalGrid(nr=nr, rmin=0., rmax=xmax, bc_rmax='reflective',
                                 nz=nz, zmin=zmin, zmax=zmax, bc_zmin='open', bc_zmax='open',
                                 n_azimuthal_modes=2,
                                 moving_window_velocity = moving_window_velocity,
                                 warpx_max_grid_size=32)

smoother = picmi.BinomialSmoother(n_pass=2, compensation=True)
solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.,
                                     method=EMsolver_method,
                                     source_smoother=smoother)

# Diagnostics
# -----------

field_diag = picmi.FieldDiagnostic(grid = grid,
                                    period = 100,
                                    warpx_plot_raw_fields = 1,
                                    warpx_plot_raw_fields_guards = 1,
                                    warpx_plot_finepatch = 1,
                                    warpx_plot_crsepatch = 1)

part_diag = picmi.ParticleDiagnostic(period = 100,
                                      species = [electrons])

# Simulation setup
# -----------------

sim = picmi.Simulation(solver=solver, verbose = 1)

sim.add_species(electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=number_per_cell_each_dim))

sim.add_laser(laser, injection_method=laser_antenna)

sim.add_diagnostic(field_diag1)
sim.add_diagnostic(part_diag1)


# Picmi input script
# ==================

run_python_simulation = True

if run_python_simulation:
    # `sim.step` will run the code, controlling it from Python
    sim.step(max_steps)
else:
    # `write_inputs` will create an input file that can be used to run
    # with the compiled version.
    sim.set_max_step(max_steps)
    sim.write_input_file(file_name='input_script')
