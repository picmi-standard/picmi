import PICMI_Base
codename = 'dummy_code'

picmi_classes = [ cls for cls in PICMI_Base.__dir__() \
                      if cls.startswith('PICMI_') ]

for cls in picmi_classes:

    def_class_code = \
"""
class %s(PICMI_Base.%s):
    def init(self):
        pass
""" %(cls[len('PICMI_'):], cls)
    exec( def_class_code, globals(), locals() )
