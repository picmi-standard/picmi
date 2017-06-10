"""Physical constants
SI units

Values from 2014CODATA

http:physics.nist.govcuuConstantsbibliography.html

"""
from .PICMI_MathConstants import pi

CODATA = '2014'

clight = 2.99792458e+8  # Speed of light in vacuum (exact) [m/s]
amu = 1.660539040e-27  # Atomic Mass Unit [kg]
echarge = 1.6021766208e-19  # Proton charge [C]
emass = 9.10938356e-31  # Electron mass [kg]
boltzmann = 1.38064852e-23  # Boltzmann's constant [J/K]
avogadro = 6.022140857e23  # Avogadro's Number [1]
planck = 6.626070040e-34  # Planck's constant [Js]

# --- Magnetic constant, permeability of free space = 4*pi*1.e-7 [N/A**2]
mu0 = 4.*pi*1.e-7

# --- Conversion factor from joules to eV, J per eV
jperev = echarge

# --- Permitivity of free space
eps0 = 1./(mu0*clight*clight)

