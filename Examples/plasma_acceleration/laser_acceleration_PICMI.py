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

# This should be the only line that needs to be changed
from fbpic import picmi
# from warp import picmi

# Define laser
laser = picmi.GaussianLaser(
    a0=laser_a0,
    wavelength=0.8e-6,  # The wavelength of the laser (in meters)
    waist=laser_waist,  # The waist of the laser (in meters)
    duration=laser_duration,  # The duration of the laser (in seconds)
    polarization_angle = np.pi/2,  # The main polarization vector
    focal_position=(0., 0., 100.e-6),  # Focal position (m)
    centroid_position=(0., 0., 0.),  # Position of the laser centroid at t=0 (m)
    propagation_direction=(0., 0., 1.) )

# Define electron beam
beam_dist = picmi.GaussianBunchDistribution(
                n_physical_particles=1.e9,
                gamma=1000,
                x_rms=2.e-6, y_rms=2.e-6, z_rms=1.e-6,
                x_emittance=1.e-6, y_emittance=1.e-6, z_emittance=1.e-6,
                z_focus=100.e-6, t_focus=10.e-15 )
                # Initializes the beam so that it will focus at z_focus,
                # at time t_focus (should we use Twiss parameters instead?)
beam = picmi.Species(
            particle_type=picmi.Electron,
            initial_distribution=beam_dist )

# Define plasma
plasma_dist = picmi.DistributionFromParsedExpression(
                density_expression="1.e23*tanh((z - 20.e-6)/100.e-6)" )
                # Plasma density in m-3
plasma = picmi.MultiSpecies(
                particle_types=[ picmi.Electron, picmi.Helium, picmi.Argon ],
                proportions=[ 1., 0.2, 0.8 ],
                initial_distribution=plasma_dist )

"""
Numerics part - can be in separate file
"""
import numpy as np
from scipy.constants import c
import picmi

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
if picmi.code == 'fbpic':
    grid = picmi.CylindricalGrid(
        nr=nx, rmin=0., rmax=xmax, bc_rmax='reflective',
        nz=nz, zmin=zmin, zmax=zmax, bc_zmin='open', bc_zmax='open',
        n_azimuthal_modes=2,
        moving_window_velocity=v_window )
elif picmi.code in ['warp', 'warpx']:
    grid_specific_arguments = {}
    if  picmi.code == 'warpx':
        grid_specific_arguments = {'max_grid_size':32, 'max_level':0}
    grid = picmi.Cartesian3DGrid(
        nx=nx, xmin=xmin, xmax=xmax, bc_xmin='periodic', bc_xmax='periodic',
        ny=ny, ymin=ymin, ymax=ymax, bc_ymin='periodic', bc_ymax='periodic',
        nz=nz, zmin=zmin, zmax=zmax, bc_zmin='open', bc_zmax='open',
        moving_window_velocity=v_window, **grid_specific_arguments )

# Setup the electromagnetic solver
solver = picmi.ElectromagneticSolver( grid=grid, cfl=1.0 )

# Initialize the simulation object
if picmi.code == 'warp':
    sim_specific_arguments = {
        'current_deposition_algo':'Esirkepov',
        'charge_deposition_algo':'Direct',
        'field_gathering_algo':'Energy_conserving',
        'particle_pusher_algo':'Boris' }
else:
    sim_specific_arguments = {
        'current_correction':'curl-free'}

sim = picmi.Simulation(
        solver=solver,
        **sim_specific_arguments,
        dt=None )  # Takes dt from the solver

# Inject the laser through an antenna
antenna = picmi.LaserAntenna(
                position=(0, 0, 9.e-6),
                normal_vector=(0, 0, 1.) )
sim.add_laser( laser, injection_method=antenna )

# Add the plasma: continuously injected by the moving window
if picmi.code == 'fbpic':
    n_macroparticle_per_cell={'r':2, 'z':2, 'theta':4}
else:
    n_macroparticle_per_cell={'x':2, 'y':2, 'z':2}

plasma_layout = picmi.EvenlySpacedLayout(
                    grid=grid,
                    n_macroparticle_per_cell=n_macroparticle_per_cell,
                    continuous_injection=True )
sim.add_species( species=plasma, layout=plasma_layout )
# For Python-driven codes: macroparticles are created at this point

# Add the beam
beam_layout = picmi.RandomDraw(
                n_macroparticles=10**5,
                seed=0 )
sim.add_species( species=beam, layout=beam_layout )

"""
picmi input script
"""
import numpy as np
from pywarpx import picmi

if picmi.code == 'warpx':
    sim.set_max_step( 1000 )
    sim.write_input_file( file_name='warpx_input_script' )
else:
    sim.step( 1000 )
