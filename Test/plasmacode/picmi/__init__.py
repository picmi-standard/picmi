import picmistandard
codename = 'plasmacode'
picmistandard.register_codename(codename)

class constants:
    # --- Put the constants in their own namespace
    c = 299792458.
    ep0 = 8.8541878128e-12
    mu0 = 1.25663706212e-06
    q_e = 1.602176634e-19
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    hbar = 1.054571817e-34
    kb = 1.380649e-23

picmistandard.register_constants(constants)

picmi_classes = [ cls for cls in picmistandard.__dir__() \
                      if cls.startswith('PICMI_') ]

for cls in picmi_classes:
    # --- Skip classes that need special handling
    if cls in ['PICMI_CylindricalGrid', 'PICMI_Simulation']:
        continue

    def_class_code = \
"""
class %s(picmistandard.%s):
    pass
""" %(cls[len('PICMI_'):], cls)
    exec( def_class_code, globals(), locals() )

picmistandard.PICMI_MultiSpecies.Species_class = Species

class CylindricalGrid(picmistandard.PICMI_CylindricalGrid):
    def init(self, kw):
        self.mode_phase = kw.pop('plasmacode_mode_phase', 0.)

class Simulation(picmistandard.PICMI_Simulation):
    def init(self, kw):
        pass
    def write_input_file(self, file_name):
        pass
    def step(self, nsteps=1):
        pass
