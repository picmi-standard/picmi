# 2D version of the LWFA example, provided for quick debugging
# This should be the only line that needs to be changed for different codes
# e.g. `from pywarpx import picmi`
#      `from fbpic import picmi`
#      `from warp import picmi`
from fbpic import picmi

# Create alias fror constants
cst = picmi.constants

# Run parameters - can be in separate file
# ========================================

# Physics parameters
# ------------------

# --- laser
import math
laser_polarization   = math.pi/2 # Polarization angle (in rad)
laser_a0             = 4.        # Normalized potential vector
laser_wavelength     = 8e-07     # Wavelength of the laser (in meters)
laser_waist          = 5e-06     # Waist of the laser (in meters)
laser_duration       = 15e-15    # Duration of the laser (in seconds)
laser_injection_loc  = 9.e-6     # Position of injection (in meters, along z)
laser_focal_distance = 100.e-6   # Focal distance from the injection (in meters)
laser_t_peak         = 30.e-15   # The time at which the laser reaches its peak
                                 # at the antenna injection location (in seconds)
# --- plasma
plasma_density_expression = "1.e23*(1+tanh((y - 20.e-6)/10.e-6))/2."
plasma_min     = [-20.e-6,  0.0e-6]
plasma_max     = [ 20.e-6,  1.e-3]

# --- electron bunch
# bunch_physical_particles  = 1.e8
# in 2D or even 1D situations where the bunch is 
# infinite in the ignorable direction, then the bunch
# density is the quantity that makes more sense
# For the example where sigma is 1 micron in all 3 directions
# and N is 10^8, this corresponds to a bunch density of
# 10^8/(1e-6**3(2*pi)**3/2) = 6.35e24
bunch_density = 6.35e24
bunch_rms_size            = [1.e-6, 1.e-6]
bunch_rms_velocity        = [0.,10.*cst.c]
bunch_centroid_position   = [0.,-22.e-6]
bunch_centroid_velocity   = [0.,1000.*cst.c]

# Numerics parameters
# -------------------

# --- Nb time steps
max_steps = 1000

# --- grid
nx = 64

ny = 480
xmin = 1.5*plasma_min[0]
xmax = 1.5*plasma_max[0]
ymin = -38.e-6
ymax = 10.e-6
moving_window_velocity = [0., cst.c]
n_macroparticle_per_cell = [4, 4]

# --- geometry and solver
em_solver_method = 'CKC'  # Cole-Karkkainen-Cowan stencil
geometry = '2D'

# Physics part - can be in separate file
# ======================================

# Physics components
# ------------------

# --- laser
laser = picmi.GaussianLaser(
    wavelength             = laser_wavelength,
    waist                  = laser_waist,
    duration               = laser_duration,
    focal_position         = [0., laser_focal_distance + laser_injection_loc],
    centroid_position      = [0., laser_injection_loc - cst.c*laser_t_peak],
    polarization_direction = [math.sin(laser_polarization), 0., math.cos(laser_polarization),],
    propagation_direction  = [0,1],
    a0                     = laser_a0)

# --- plasma
plasma_dist = picmi.AnalyticDistribution(
                density_expression = plasma_density_expression,
                lower_bound        = plasma_min,
                upper_bound        = plasma_max,
                fill_in            = True)
plasma = picmi.MultiSpecies(
                particle_types = ['He', 'Ar', 'electron'],
                names          = ['He+', 'Argon', 'e-'],
                charge_states  = [1, 5, None],
                proportions    = [0.2, 0.8, 0.2 + 5*0.8],
                initial_distribution=plasma_dist)

# --- electron bunch
beam_dist = picmi.GaussianBunchDistribution(
#            n_physical_particles = bunch_physical_particles,
            density = bunch_density
            rms_bunch_size       = bunch_rms_size,
            rms_velocity         = bunch_rms_velocity,
            centroid_position    = bunch_centroid_position,
            centroid_velocity    = bunch_centroid_velocity )
beam = picmi.Species( particle_type        = 'electron',
                      name                 = 'beam',
                      initial_distribution = beam_dist)

# Numerics components
# -------------------

grid = picmi.Cartesian2DGrid(
    number_of_cells           = [nx, ny],
    lower_bound               = [xmin, ymin],
    upper_bound               = [xmax, ymax],
    lower_boundary_conditions = ['periodic', 'open'],
    upper_boundary_conditions = ['periodic', 'open'],
    moving_window_velocity    = moving_window_velocity,
    )


smoother = picmi.BinomialSmoother( n_pass       = 1,
                                   compensation = True )
solver = picmi.ElectromagneticSolver( grid            = grid,
                                      cfl             = 1.,
                                      method          = em_solver_method,
                                      source_smoother = smoother)

# Diagnostics
# -----------
field_diag = picmi.FieldDiagnostic(name = 'diag1',
                                   grid = grid,
                                   period = 100,)
part_diag = picmi.ParticleDiagnostic(name = 'diag1',
                                     period = 100,
                                     species = [beam])

# Simulation setup
# -----------------

# Initialize the simulation object
# Note that the time step size is obtained from the solver
sim = picmi.Simulation(solver = solver, verbose = 1)

# Inject the laser through an antenna
antenna = picmi.LaserAntenna(
                position = (0, 0, 9.e-6),
                normal_vector = (0, 0, 1.))
sim.add_laser(laser, injection_method = antenna)

# Add the plasma: continuously injected by the moving window
plasma_layout = picmi.GriddedLayout(
                    grid = grid,
                    n_macroparticle_per_cell = n_macroparticle_per_cell)
sim.add_species(species=plasma, layout=plasma_layout)

# Add the beam
beam_layout = picmi.PseudoRandomLayout(
                n_macroparticles = 10**5,
                seed = 0)
initialize_self_field = True
if picmi.codename == 'warpx':
    initialize_self_field = False
sim.add_species(species=beam, layout=beam_layout,
                initialize_self_field=initialize_self_field)

# Add the diagnostics
sim.add_diagnostic(field_diag)
sim.add_diagnostic(part_diag)

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
