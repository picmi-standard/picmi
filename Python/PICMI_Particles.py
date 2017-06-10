from .CODATA2014 import *

class Particle(object):
    def __init__(self, mass=None, charge=None, symbol=None, name=None):
        self.symbol = symbol
        self.name = name
        if mass is not None:
            self.mass = mass
            self.M = self.mass
        if charge is not None:
            self.charge = charge
            self.Q = self.charge


class Atom(Particle):
    def __init__(self, symbol, A, Z, group, period, name=None):
        Particle.__init__(self, mass=A*amu, symbol=symbol, name=name)
        self.A = A
        self.Z = Z
        self.group = group
        self.period = period


class Molecule(Particle):
    def __init__(self, symbol, composition, name=None):
        self.composition = composition
        mass = 0.
        self.A = 0
        self.Z = 0
        for c in composition:
            mass += c.mass
            self.A += c.A
            self.Z += c.Z
        Particle.__init__(self, mass=mass, symbol=symbol, name=name)


Hydrogen = Atom(A=1.00794, symbol='H', Z=1, group=1, period=1, name='Hydrogen')
Deuterium = Atom(A=2.01410178, symbol='D', Z=1, group=1, period=1, name='Deuterium')
Tritium = Atom(A=3.0160492, symbol='T', Z=1, group=1, period=1, name='Tritium')
Helium = Atom(A=4.002602, symbol='He', Z=2, group=18, period=1, name='Helium')
Lithium = Atom(A=6.941, symbol='Li', Z=3, group=1, period=2, name='Lithium')
Lithium6 = Atom(A=6.01512279516, symbol='Li6', Z=3, group=1, period=2, name='Lithium6')
Lithium7 = Atom(A=7.016004558, symbol='Li7', Z=3, group=1, period=2, name='Lithium7')
Beryllium = Atom(A=9.012182, symbol='Be', Z=4, group=2, period=2, name='Beryllium')
Boron = Atom(A=10.811, symbol='B', Z=5, group=13, period=2, name='Boron')
Carbon = Atom(A=12.0107, symbol='C', Z=6, group=14, period=2, name='Carbon')
Nitrogen = Atom(A=14.0067, symbol='N', Z=7, group=15, period=2, name='Nitrogen')
Oxygen = Atom(A=15.9994, symbol='O', Z=8, group=16, period=2, name='Oxygen')
Fluorine = Atom(A=18.9984032, symbol='F', Z=9, group=17, period=2, name='Fluorine')
Neon = Atom(A=20.1797, symbol='Ne', Z=10, group=18, period=2, name='Neon')
Sodium = Atom(A=22.98977, symbol='Na', Z=11, group=1, period=3, name='Sodium')
Magnesium = Atom(A=24.305, symbol='Mg', Z=12, group=2, period=3, name='Magnesium')
Aluminium = Atom(A=26.981538, symbol='Al', Z=13, group=13, period=3, name='Aluminium')
Silicon = Atom(A=28.0855, symbol='Si', Z=14, group=14, period=3, name='Silicon')
Phosphorus = Atom(A=30.973761, symbol='P', Z=15, group=15, period=3, name='Phosphorus')
Sulfur = Atom(A=32.065, symbol='S', Z=16, group=16, period=3, name='Sulfur')
Chlorine = Atom(A=35.453, symbol='Cl', Z=17, group=17, period=3, name='Chlorine')
Argon = Atom(A=39.948, symbol='Ar', Z=18, group=18, period=3, name='Argon')
Potassium = Atom(A=39.0983, symbol='K', Z=19, group=1, period=4, name='Potassium')
Calcium = Atom(A=40.078, symbol='Ca', Z=20, group=2, period=4, name='Calcium')
Scandium = Atom(A=44.95591, symbol='Sc', Z=21, group=3, period=4, name='Scandium')
Titanium = Atom(A=47.867, symbol='Ti', Z=22, group=4, period=4, name='Titanium')
Vanadium = Atom(A=50.9415, symbol='V', Z=23, group=5, period=4, name='Vanadium')
Chromium = Atom(A=51.9961, symbol='Cr', Z=24, group=6, period=4, name='Chromium')
Manganese = Atom(A=54.938049, symbol='Mn', Z=25, group=7, period=4, name='Manganese')
Iron = Atom(A=55.845, symbol='Fe', Z=26, group=8, period=4, name='Iron')
Cobalt = Atom(A=58.9332, symbol='Co', Z=27, group=9, period=4, name='Cobalt')
Nickel = Atom(A=58.6934, symbol='Ni', Z=28, group=10, period=4, name='Nickel')
Copper = Atom(A=63.546, symbol='Cu', Z=29, group=11, period=4, name='Copper')
Zinc = Atom(A=65.409, symbol='Zn', Z=30, group=12, period=4, name='Zinc')
Gallium = Atom(A=69.723, symbol='Ga', Z=31, group=13, period=4, name='Gallium')
Germanium = Atom(A=72.64, symbol='Ge', Z=32, group=14, period=4, name='Germanium')
Arsenic = Atom(A=74.9216, symbol='As', Z=33, group=15, period=4, name='Arsenic')
Selenium = Atom(A=78.96, symbol='Se', Z=34, group=16, period=4, name='Selenium')
Bromine = Atom(A=79.904, symbol='Br', Z=35, group=17, period=4, name='Bromine')
Krypton = Atom(A=83.798, symbol='Kr', Z=36, group=18, period=4, name='Krypton')
Rubidium = Atom(A=85.4678, symbol='Rb', Z=37, group=1, period=5, name='Rubidium')
Strontium = Atom(A=87.62, symbol='Sr', Z=38, group=2, period=5, name='Strontium')
Yttrium = Atom(A=88.90585, symbol='Y', Z=39, group=3, period=5, name='Yttrium')
Zirconium = Atom(A=91.224, symbol='Zr', Z=40, group=4, period=5, name='Zirconium')
Niobium = Atom(A=92.90638, symbol='Nb', Z=41, group=5, period=5, name='Niobium')
Molybdenum = Atom(A=95.94, symbol='Mo', Z=42, group=6, period=5, name='Molybdenum')
Technetium = Atom(A=98.0, symbol='Tc', Z=43, group=7, period=5, name='Technetium')
Ruthenium = Atom(A=101.07, symbol='Ru', Z=44, group=8, period=5, name='Ruthenium')
Rhodium = Atom(A=102.9055, symbol='Rh', Z=45, group=9, period=5, name='Rhodium')
Palladium = Atom(A=106.42, symbol='Pd', Z=46, group=10, period=5, name='Palladium')
Silver = Atom(A=107.8682, symbol='Ag', Z=47, group=11, period=5, name='Silver')
Cadmium = Atom(A=112.411, symbol='Cd', Z=48, group=12, period=5, name='Cadmium')
Indium = Atom(A=114.818, symbol='In', Z=49, group=13, period=5, name='Indium')
Tin = Atom(A=118.71, symbol='Sn', Z=50, group=14, period=5, name='Tin')
Antimony = Atom(A=121.76, symbol='Sb', Z=51, group=15, period=5, name='Antimony')
Tellurium = Atom(A=127.6, symbol='Te', Z=52, group=16, period=5, name='Tellurium')
Iodine = Atom(A=126.90447, symbol='I', Z=53, group=17, period=5, name='Iodine')
Xenon = Atom(A=131.293, symbol='Xe', Z=54, group=18, period=5, name='Xenon')
Caesium = Atom(A=132.90545, symbol='Cs', Z=55, group=1, period=6, name='Caesium')
Cesium = Caesium
Barium = Atom(A=137.327, symbol='Ba', Z=56, group=2, period=6, name='Barium')
Lutetium = Atom(A=174.967, symbol='Lu', Z=71, group=3, period=6, name='Lutetium')
Hafnium = Atom(A=178.49, symbol='Hf', Z=72, group=4, period=6, name='Hafnium')
Tantalum = Atom(A=180.9479, symbol='Ta', Z=73, group=5, period=6, name='Tantalum')
Tungsten = Atom(A=183.84, symbol='W', Z=74, group=6, period=6, name='Tungsten')
Rhenium = Atom(A=186.207, symbol='Re', Z=75, group=7, period=6, name='Rhenium')
Osmium = Atom(A=190.23, symbol='Os', Z=76, group=8, period=6, name='Osmium')
Iridium = Atom(A=192.217, symbol='Ir', Z=77, group=9, period=6, name='Iridium')
Platinum = Atom(A=195.078, symbol='Pt', Z=78, group=10, period=6, name='Platinum')
Gold = Atom(A=196.96655, symbol='Au', Z=79, group=11, period=6, name='Gold')
Mercury = Atom(A=200.59, symbol='Hg', Z=80, group=12, period=6, name='Mercury')
Thallium = Atom(A=204.3833, symbol='Tl', Z=81, group=13, period=6, name='Thallium')
Lead = Atom(A=207.2, symbol='Pb', Z=82, group=14, period=6, name='Lead')
Bismuth = Atom(A=208.98038, symbol='Bi', Z=83, group=15, period=6, name='Bismuth')
Polonium = Atom(A=210.0, symbol='Po', Z=84, group=16, period=6, name='Polonium')
Astatine = Atom(A=210.0, symbol='At', Z=85, group=17, period=6, name='Astatine')
Radon = Atom(A=220.0, symbol='Rn', Z=86, group=18, period=6, name='Radon')
Francium = Atom(A=223.0, symbol='Fr', Z=87, group=1, period=7, name='Francium')
Radium = Atom(A=226.0, symbol='Ra', Z=88, group=2, period=7, name='Radium')
Actinium = Atom(A=227.0, symbol='Ac', Z=89, group=3, period=7, name='Actinium')
Thorium = Atom(A=232.0381, symbol='Th', Z=90, group=102, period=7, name='Thorium')
Protactinium = Atom(A=231.0359, symbol='Pa', Z=91, group=102, period=7, name='Protactinium')
Uranium = Atom(A=238.02891, symbol='U', Z=92, group=102, period=7, name='Uranium')
Neptunium = Atom(A=237.0, symbol='Np', Z=93, group=102, period=7, name='Neptunium')
Plutonium = Atom(A=244.0, symbol='Pu', Z=94, group=102, period=7, name='Plutonium')
Americium = Atom(A=243.0, symbol='Am', Z=95, group=102, period=7, name='Americium')
Curium = Atom(A=247.0, symbol='Cm', Z=96, group=102, period=7, name='Curium')
Berkelium = Atom(A=247.0, symbol='Bk', Z=97, group=102, period=7, name='Berkelium')
Californium = Atom(A=251.0, symbol='Cf', Z=98, group=102, period=7, name='Californium')
Einsteinium = Atom(A=252.0, symbol='Es', Z=99, group=102, period=7, name='Einsteinium')
Fermium = Atom(A=257.0, symbol='Fm', Z=100, group=102, period=7, name='Fermium')
Mendelevium = Atom(A=258.0, symbol='Md', Z=101, group=102, period=7, name='Mendelevium')
Nobelium = Atom(A=259.0, symbol='No', Z=102, group=102, period=7, name='Nobelium')
Lawrencium = Atom(A=262.0, symbol='Lr', Z=103, group=102, period=7, name='Lawrencium')
Rutherfordium = Atom(A=261.0, symbol='Rf', Z=104, group=4, period=7, name='Rutherfordium')
Lawrencium = Atom(A=262.0, symbol='Lr', Z=103, group=3, period=7, name='Lawrencium')
Dubnium = Atom(A=262.0, symbol='Db', Z=105, group=5, period=7, name='Dubnium')
Bohrium = Atom(A=264.0, symbol='Bh', Z=107, group=7, period=7, name='Bohrium')
Seaborgium = Atom(A=266.0, symbol='Sg', Z=106, group=6, period=7, name='Seaborgium')
Meitnerium = Atom(A=268.0, symbol='Mt', Z=109, group=9, period=7, name='Meitnerium')
Darmstadtium = Atom(A=271.0, symbol='Ds', Z=110, group=10, period=7, name='Darmstadtium')
Roentgenium = Atom(A=272.0, symbol='Rg', Z=111, group=11, period=7, name='Roentgenium')
Hassium = Atom(A=277.0, symbol='Hs', Z=108, group=8, period=7, name='Hassium')
Meitnerium = Atom(A=278.0, symbol='Mt', Z=109, group=9, period=7, name='Meitnerium')
Darmstadtium = Atom(A=281.0, symbol='Ds', Z=110, group=10, period=7, name='Darmstadtium')
Roentgenium = Atom(A=282.0, symbol='Rg', Z=111, group=11, period=7, name='Roentgenium')
Copernicium = Atom(A=285.0, symbol='Cn', Z=112, group=12, period=7, name='Copernicium')
Nihonium = Atom(A=286.0, symbol='Nh', Z=113, group=13, period=7, name='Nihonium')
Flerovium = Atom(A=289.0, symbol='Fl', Z=114, group=14, period=7, name='Flerovium')
Moscovium = Atom(A=290.0, symbol='Mc', Z=115, group=15, period=7, name='Moscovium')
Livermorium = Atom(A=293.0, symbol='Lv', Z=116, group=16, period=7, name='Livermorium')
Tennessine = Atom(A=294.0, symbol='Ts', Z=117, group=17, period=7, name='Tennessine')
Oganesson = Atom(A=294.0, symbol='Og', Z=118, group=18, period=7, name='Oganesson')

Electron = Particle(charge=-echarge, mass=emass, symbol='e-', name='Electron')
Positron = Particle(charge=+echarge, mass=emass, symbol='e+', name='Positron')
Muon = Particle(charge=-echarge, mass=1.883531475e-28, symbol='mu-', name='Muon')
Antimuon = Particle(charge=+echarge, mass=1.883531475e-28, symbol='mu+', name='Antimuon')

Proton = Particle(mass=1.6726231e-27, charge=echarge, symbol='P', name='Proton')
AntiProton = Particle(mass=1.6726231e-27, charge=-echarge, symbol='P', name='Antiproton')
Neutron = Particle(mass=1.6749286e-27, charge=0, symbol='N', name='Neutron')

Dihydrogen = Molecule(composition=[Hydrogen, Hydrogen], symbol='H2', name='Dihydrogen')
Dideuterium = Molecule(composition=[Deuterium, Deuterium], symbol='D2', name='Dideuterium')
Dinitrogen = Molecule(composition=[Nitrogen, Nitrogen], symbol='N2', name='Dinitrogen')
Dioxygen = Molecule(composition=[Oxygen, Oxygen], symbol='O2', name='Dioxygen')
Carbon_Monoxide = Molecule(composition=[Carbon, Oxygen], symbol='CO', name='Carbon_Monoxide')
Carbon_Dioxide = Molecule(composition=[Carbon, Oxygen, Oxygen], symbol='CO2', name='Carbon_Dioxide')
Water = Molecule(composition=[Hydrogen, Hydrogen, Oxygen], symbol='H2O', name='Water')

Photon = Particle(charge=0., mass=0., symbol='gnu', name='Photon')

Hydrogen.ionization_levels = [ 13.6]

Helium.ionization_levels = [ 24.58,
                             54.4]

Lithium.ionization_levels = [ 5.392,
                              75.638,
                              122.451]

Boron.ionization_levels = [ 8.297,
                            25.154,
                            37.929,
                            259.367,
                            340.216]

Carbon.ionization_levels = [ 11.261,
                             24.383,
                             47.887,
                             64.492,
                             392.081,
                             489.979]

Nitrogen.ionization_levels = [ 14.534,
                               29.601,
                               47.448,
                               77.472,
                               97.888,
                               552.057,
                               667.029]

Neon.ionization_levels = [ 21.564,
                           40.962,
                           63.45,
                           97.11,
                           126.21,
                           157.93,
                           207.27,
                           239.09,
                           1195.79,
                           1362.164]

Aluminium.ionization_levels = [ 5.968,
                                18.796,
                                28.399,
                                119.78,
                                153.561,
                                190.156,
                                241.34,
                                284.163,
                                329.564,
                                398.057,
                                441.232,
                                2082.379,
                                2300.16]

Silicon.ionization_levels = [ 7.264,
                              16.95,
                              34.27,
                              46.65,
                              159.8,
                              210.5,
                              261.3,
                              312.0,
                              364.0,
                              415.1,
                              503.7,
                              552.2,
                              2324.0,
                              2569.]

Argon.ionization_levels = [ 15.759,
                            27.6,
                            40.7,
                            59.81,
                            75.02,
                            91.007,
                            124.319,
                            143.456,
                            422.44,
                            478.68,
                            538.95,
                            618.24,
                            686.09,
                            755.73,
                            854.75,
                            918.0,
                            4120.778,
                            4426.114]

Copper.ionization_levels = [ 7.70,
                             20.234,
                             36.739,
                             57.213,
                             79.577,
                             102.313,
                             138.48,
                             165.354,
                             198.425,
                             231.492,
                             264.566,
                             367.913,
                             399.951,
                             434.055,
                             482.627,
                             518.798,
                             554.970,
                             631.446,
                             668.672,
                             1691.78,
                             1799.261,
                             1916.0,
                             2060.0,
                             2182.0,
                             2308.0,
                             2478.0,
                             2587.5,
                             11062.38,
                             11567.617]

Krypton.ionization_levels = [ 13.999,
                              24.359,
                              36.95,
                              52.5,
                              64.7,
                              78.5,
                              111.0,
                              125.8,
                              230.85,
                              268.2,
                              308.0,
                              350.0,
                              391.0,
                              447.0,
                              492.0,
                              541.0,
                              592.0,
                              641.0,
                              786.0,
                              833.0,
                              884.0,
                              937.0,
                              998.0,
                              1051.0,
                              1151.0,
                              1205.3,
                              2928.0,
                              3070.0,
                              3227.0,
                              3381.]

# The Xenon data is taken from NIST (http://physics.nist.gov/)
Xenon.ionization_levels = [ 12.1298431,
                            20.975,
                            31.05,
                            42.20,
                            54.1,
                            66.703,
                            91.6,
                            105.9778,
                            179.84,
                            202.0,
                            229.02,
                            255.0,
                            280.9,
                            314.1,
                            342.9,
                            374.1,
                            403.9,
                            433.9,
                            549.0,
                            582.0,
                            616.0,
                            650.0,
                            700.0,
                            736.0,
                            818.0,
                            857.0,
                            1493.0,
                            1571.0,
                            1653.0,
                            1742.0,
                            1826.0,
                            1919.0,
                            2023.0,
                            2113.0,
                            2209.0,
                            2300.0,
                            2556.0,
                            2637.0,
                            2726.0,
                            2811.0,
                            2975.0,
                            3068.0,
                            3243.0,
                            3333.8,
                            7660.0,
                            7889.0,
                            8144.0,
                            8382.0,
                            8971.0,
                            9243.0,
                            9581.0,
                            9810.0,
                            40272.0,
                            41300.0]

Cesium.ionization_levels = [ 3.894]

