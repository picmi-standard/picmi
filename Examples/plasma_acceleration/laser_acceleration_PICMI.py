"""
Run parameters - can be in separate file
"""
# Laser parameters
laser_waist = 5.e-6  # The waist of the laser (in meters)
laser_duration = 15.e-15  # The duration of the laser (in seconds)
laser_a0 = 4. # Amplitude of the normalized vector potential

"""
Physics part - can be in separate file
"""
import numpy as np

# This should be the only line that needs to be changed for different codes
from plasma_code import picmi

# Define the laser profile laser
laser = picmi.GaussianLaser(
            a0 = laser_a0,
            wavelength = 0.8e-6,  # The wavelength of the laser (in meters)
            waist = laser_waist,  # The waist of the laser (in meters)
            duration = laser_duration,  # The duration of the laser (in seconds)
            polarization_angle = np.pi/2,  # The main polarization vector
            focal_position = [0., 0., 100.e-6],  # Focal position (m)
            centroid_position = [0., 0., 0.],  # Position of the laser centroid at t=0 (m)
            propagation_direction = [0., 0., 1.])

# Define the electron beam
beam_dist = picmi.GaussianBunchDistribution(
                n_physical_particles = 1.e9,
                rms_bunch_size = [2.e-6, 2.e-6, 1.e-6],
                rms_velocity = [0.,0.,10.],
                centroid_position = [0.,0.,100.e-6],
                centroid_velocity = [0.,0.,1000.],
                velocity_divergence = [0.,0.,0.])

beam = picmi.Species(
            particle_type = 'electron',
            initial_distribution = beam_dist)

# Define plasma
plasma_dist = picmi.AnalyticDistribution(
                density_expression = "1.e23*tanh((z - 20.e-6)/100.e-6)",
                fill_in=True)

plasma = picmi.MultiSpecies(
                particle_types = ['He', 'Ar', 'electron'],
                names = ['He+', 'Argon', 'e-'],
                charge_states = [1, 5, None],
                proportions = [0.2, 0.8, 0.2 + 5*0.8],
                initial_distribution=plasma_dist)

# Individual species in a `MultiSpecies` can be addressed either
# with their index (using Python indexing conventions) or with their name
# (if the user provided a name)
# Set the ionization for the species number 1 (Argon)
# and place the created electrons into the species number 2 (electron)
plasma['Argon'].activate_ionization(model="ADK", target_species=plasma['e-'])

"""
Numerics part - can be in separate file
"""
import numpy as np
from scipy.constants import c
from plasma_code import picmi

# Define the grid
nx = 64
ny = 64
nz = 480
xmin = -30.e-6
xmax = +30.e-6
ymin = -30.e-6
ymax = +30.e-6
zmin = -56.e-6
zmax = +12.e-6
v_window = (0., 0., c)

# Setup the grid ; this may be code dependent
# Note that code specific arguments use the code name as a prefix.
if picmi.codename == 'plasma_code':
    nr = nx/2.
    grid = picmi.CylindricalGrid(
        nr=nr, rmin=0., rmax=xmax, bc_rmax='reflective',
        nz=nz, zmin=zmin, zmax=zmax, bc_zmin='open', bc_zmax='open',
        n_azimuthal_modes=2,
        moving_window_velocity=v_window,
        plasma_code_mode_phase = 0.)
elif picmi.codename in ['plasma_code_3D', 'AMRplasma_code']:
    grid = picmi.Cartesian3DGrid(
        nx=nx, xmin=xmin, xmax=xmax, bc_xmin='periodic', bc_xmax='periodic',
        ny=ny, ymin=ymin, ymax=ymax, bc_ymin='periodic', bc_ymax='periodic',
        nz=nz, zmin=zmin, zmax=zmax, bc_zmin='open', bc_zmax='open',
        moving_window_velocity=v_window,
        AMRplasma_code_max_grid_size=32, AMRplasma_code_max_level=0)

# Setup the electromagnetic solver
smoother = picmi.BinomialSmoother(n_pass=2, compensation=True)
solver = picmi.ElectromagneticSolver(grid=grid, cfl=1.0,
                                     source_smoother=smoother)

# Initialize the simulation object
# Note that the time step size is obtained from the solver
sim = picmi.Simulation(
        solver=solver)

# Inject the laser through an antenna
antenna = picmi.LaserAntenna(
                position = (0, 0, 9.e-6),
                normal_vector = (0, 0, 1.))
sim.add_laser(laser, injection_method = antenna)

# Add the plasma: continuously injected by the moving window
if picmi.codename == 'plasma_code':
    n_macroparticle_per_cell = {'r':2, 'z':2, 'theta':4}
else:
    n_macroparticle_per_cell = {'x':2, 'y':2, 'z':2}

plasma_layout = picmi.GriddedLayout(
                    grid = grid,
                    n_macroparticle_per_cell = n_macroparticle_per_cell)
sim.add_species(species=plasma, layout=plasma_layout)
# For Python-driven codes: macroparticles are created at this point

# Add the beam
beam_layout = picmi.PseudoRandomLayout(
                n_macroparticles = 10**5,
                seed = 0)
sim.add_species(species=beam, layout=beam_layout, initialize_self_field=True)

"""
picmi input script
"""
import numpy as np
from plasma_code import picmi

run_python_simulation = True

if run_python_simulation:
    sim.step(1000)
else:
    sim.set_max_step(1000)
    sim.write_input_file(file_name='input_script')
