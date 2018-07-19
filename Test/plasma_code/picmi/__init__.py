import PICMI_Base
code = 'plasma_code'

picmi_classes = [ cls for cls in PICMI_Base.__dir__() \
                      if cls.startswith('PICMI_') ]

for cls in picmi_classes:
    if cls == 'PICMI_Simulation':
        continue

    def_class_code = \
"""
class %s(PICMI_Base.%s):
    def init(self, **kw):
        pass
""" %(cls[len('PICMI_'):], cls)
    exec( def_class_code, globals(), locals() )

PICMI_Base.PICMI_MultiSpecies.Species_class = Species

class Simulation(PICMI_Base.PICMI_Simulation):
    def init(self, **kw):
        pass
    def write_input_file(self, file_name):
        pass
    def step(self, nsteps=1):
        pass
