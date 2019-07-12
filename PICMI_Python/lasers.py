"""Classes following the PICMI standard
These should be the base classes for Python implementation of the PICMI standard
"""
import math
import sys

from .base import _ClassWithInit

# ---------------
# Physics objects
# ---------------

class PICMI_GaussianLaser(_ClassWithInit):
    """
    Specifies a Gaussian laser distribution
      - wavelength: Laser wavelength
      - waist: Waist of the Gaussian pulse at focus [m]
      - duration: Duration of the Gaussian pulse [s]
      - focal_position=[0,0,0]: Position of the laser focus (vector) [m]
      - centroid_position=[0,0,0]: Position of the laser centroid at time 0 (vector) [m]
      - propagation_direction=[0,0,1]: Direction of propagation (unit vector) [1]
      - polarization_direction=[1,0,0]: Direction of polarization (unit vector) [1]
      - a0: Normalized vector potential at focus
            Specify either a0 or E0 (E0 takes precedence).
      - E0: Maximum amplitude of the laser field [V/m]
            Specify either a0 or E0 (E0 takes precedence).
      - zeta: Spatial chirp at focus (in the lab frame) [m.s]
      - beta: Angular dispersion at focus (in the lab frame) [rad.s]
      - phi2: Temporal chirp at focus (in the lab frame) [s^2]
    """
    def __init__(self, wavelength, waist, duration,
                 focal_position = [0., 0., 0.],
                 centroid_position = [0., 0., 0.],
                 propagation_direction = [0., 0., 1.],
                 polarization_direction = [1., 0., 0.],
                 a0 = None, 
                 E0 = None,
                 zeta = None,
                 beta = None,
                 phi2 = None,
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
        self.polarization_direction = polarization_direction
        self.a0 = a0
        self.E0 = E0
        self.zeta = zeta
        self.beta = beta
        self.phi2 = phi2

        self.handle_init(kw)


# ------------------
# Numeric Objects
# ------------------


class PICMI_LaserAntenna(_ClassWithInit):
    """
    Specifies the laser antenna injection method
      - position: Position of antenna launching the laser (vector) [m]
      - normal_vector: Vector normal to antenna plane (vector) [1]
    """
    def __init__(self, position, normal_vector, **kw):

        self.position = position
        self.normal_vector = normal_vector

        self.handle_init(kw)

