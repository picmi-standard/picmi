import PICMI_Base
codename = 'plasmacode'
PICMI_Base.register_codename(codename)

picmi_classes = [ cls for cls in PICMI_Base.__dir__() \
                      if cls.startswith('PICMI_') ]

for cls in picmi_classes:
    # --- Skip classes that need special handling
    if cls in ['PICMI_CylindricalGrid', 'PICMI_Simulation']:
        continue

    def_class_code = \
"""
class %s(PICMI_Base.%s):
    pass
""" %(cls[len('PICMI_'):], cls)
    exec( def_class_code, globals(), locals() )

PICMI_Base.PICMI_MultiSpecies.Species_class = Species

class CylindricalGrid(PICMI_Base.PICMI_CylindricalGrid):
    def init(self, kw):
        self.mode_phase = kw.pop('plasmacode_mode_phase', 0.)

class Simulation(PICMI_Base.PICMI_Simulation):
    def init(self, kw):
        pass
    def write_input_file(self, file_name):
        pass
    def step(self, nsteps=1):
        pass
