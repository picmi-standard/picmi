#!/usr/bin/env python3

# This should be the only line that needs to be changed for different codes
# e.g. `from pywarpx import picmi`
#      `from fbpic import picmi`
#      `from warp import picmi`
from plasma_code import picmi
# Create alias fror constants
cst = picmi.constants

gamma_boost = 10.                 # Lorentz factor of the boosted-frame

# Physics parameters
# ------------------

# --- laser
laser_a0              = 1.        # Normalized potential vector
laser_wavelength      = 8e-07     # Wavelength of the laser (in meters)
laser_waist           = 8e-06     # Waist of the laser (in meters)
laser_length          = 3.e-6     # Length of the laser (in meters)
laser_z0              = -6.e-6    # Initial position of the centroid (meters)
laser_focal_position  = 0.        # Laser focal position (meters)
laser_antenna_z       = -0.1e-6   # Position of the antenna (meters)

# --- plasma
plasma_density = 2.5e25
plasma_min     = [-25.e-6, -25.e-6,  5.0e-6]
plasma_max     = [ 25.e-6,  25.e-6,  3000.e-6]

# --- Beam
beam_density = 1.e26
beam_uz = 1000.
beam_xmin = -3.e-6
beam_xmax = +3.e-6
beam_ymin = -3.e-6
beam_ymax = +3.e-6
beam_zmin = -12.e-6
beam_zmax = -10.e-6

# Numerics parameters
# -------------------

# --- Nb time steps
max_steps = 500

# --- grid
nx = 100
ny = 100
nz = 800
xmin = -30.e-6
xmax = +30.e-6
ymin = -30.e-6
ymax = +30.e-6
zmin = -15.e-6
zmax = +5.e-6
moving_window_velocity = [0., 0., cst.c]
plasma_number_per_cell_each_dim = [2, 2, 2]
beam_number_per_cell_each_dim = [4, 4, 4]

# --- geometry and solver
em_solver_method = 'PSATD'
geometry = '3D'

# Note that code-specific changes can be introduced with `picmi.codename`
if picmi.codename == 'fbpic':
    geometry = 'RZ'

# Physics components
# ------------------

# --- laser
laser = picmi.GaussianLaser(
    wavelength            = laser_wavelength,
    waist                 = laser_waist,
    duration              = laser_length/cst.c,
    focal_position        = [0., 0., laser_focal_position],
    centroid_position     = [0., 0., laser_z0],
    propagation_direction = [0,0,1],
    a0                    = laser_a0)

# --- plasma
uniform_plasma = picmi.UniformDistribution(
                    density     = plasma_density,
                    lower_bound = plasma_min,
                    upper_bound = plasma_max,
                    fill_in     = True)
electrons = picmi.Species(
                    particle_type = 'electron',
                    particle_shape = 'cubic',
                    name = 'electrons',
                    initial_distribution = uniform_plasma)
ions = picmi.Species(particle_type = 'H',
                    charge_state = +1,
                    particle_shape = 'cubic',
                    name = 'ions',
                    initial_distribution = uniform_plasma)

# --- beam
beam_dist = picmi.AnalyticDistribution(
    density_expression = "beam_density* (z - beam_zmin)*(beam_zmax - z)*4/(beam_zmax - beam_zmin)**2*(1. - (sqrt(x**2 + y**2)/beam_rmax)**2)",
    beam_density = beam_density,
    beam_rmax = beam_xmax,
    beam_zmin = beam_zmin,
    beam_zmax = beam_zmax,
    lower_bound = [beam_xmin, beam_ymin, beam_zmin],
    upper_bound = [beam_xmax, beam_ymax, beam_zmax],
    directed_velocity = [0., 0., beam_uz*cst.c])
beam = picmi.Species(particle_type = 'electron',
                     particle_shape = 'cubic',
                     name = 'beam',
                     initial_distribution = beam_dist)

# Numerics components
# -------------------
if geometry == '3D':
    grid = picmi.Cartesian3DGrid(
        number_of_cells = [nx, ny, nz],
        lower_bound = [xmin, ymin, zmin],
        upper_bound = [xmax, ymax, zmax],
        lower_boundary_conditions = ['periodic', 'periodic', 'open'],
        upper_boundary_conditions = ['periodic', 'periodic', 'open'],
        moving_window_velocity = moving_window_velocity)
elif geometry == 'RZ':
    plasma_number_per_cell_each_dim[1] = 4
    beam_number_per_cell_each_dim[1] = 4
    grid = picmi.CylindricalGrid(
        number_of_cells           = [nx//2, nz],
        lower_bound               = [0., zmin],
        upper_bound               = [xmax, zmax],
        lower_boundary_conditions = [ None, 'open'],
        upper_boundary_conditions = ['reflective', 'open'],
        n_azimuthal_modes         = 2,
        moving_window_velocity    = [0, moving_window_velocity[-1]])

smoother = picmi.BinomialSmoother( n_pass = [1, 1, 1],
                                   compensation = [False, False, False] )
galilean_velocity = [0, 0, -cst.c*(1.-1./gamma_boost**2)**.5]
solver = picmi.ElectromagneticSolver(
        grid             = grid,
        method           = em_solver_method,
        cfl              = 1.,
        source_smoother  = smoother,
        galilean_velocity = galilean_velocity )

# Diagnostics
# -----------
field_diag = picmi.FieldDiagnostic(grid = grid,
                                   period = 100,
                                   write_dir = 'diags')

part_diag = picmi.ParticleDiagnostic(period = 100,
                                     species = [electrons, ions, beam],
                                     write_dir = 'diags')

field_diag_lab = picmi.LabFrameFieldDiagnostic(grid = grid,
                                               num_snapshots = 20,
                                               dt_snapshots = 0.5*(zmax - zmin)/cst.c,
                                               data_list = ["rho", "E", "B", "J"],
                                               write_dir = 'lab_diags')

part_diag_lab = picmi.LabFrameParticleDiagnostic(grid = grid,
                                                 num_snapshots = 20,
                                                 dt_snapshots = 0.5*(zmax - zmin)/cst.c,
                                                 species = [electrons, ions, beam],
                                                 write_dir = 'lab_diags')


# Simulation setup
# -----------------

sim = picmi.Simulation(solver = solver,
                       max_steps = max_steps,
                       gamma_boost = gamma_boost,
                       verbose = 1 )

sim.add_species(electrons, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=plasma_number_per_cell_each_dim))
sim.add_species(ions, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=plasma_number_per_cell_each_dim))
sim.add_species(beam, layout=picmi.GriddedLayout(grid=grid, n_macroparticle_per_cell=beam_number_per_cell_each_dim),
                                                 initialize_self_field=True)

laser_antenna = picmi.LaserAntenna(position = [0., 0., laser_antenna_z],  # This point is on the laser plane
                                   normal_vector = [0., 0., 1.])  # The plane normal direction
sim.add_laser(laser, injection_method=laser_antenna)

sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)
sim.add_diagnostic(field_diag_lab)
sim.add_diagnostic(part_diag_lab)

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
