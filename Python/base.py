"""base code for the PICMI standard
"""

codename = None

def register_codename(_codename):
    global codename
    codename = _codename

class _ClassWithInit(object):
    def handle_init(self, kw):
        codekw = {}
        for k,v in kw.items():
            if k.startswith(codename + '_'):
                codekw[k] = v

        self.init(codekw)
        if codekw:
            raise TypeError("got an unexpected keyword argument '%s'"%list(codekw)[0])

    def init(self, kw):
        # --- The implementation of this routine should use kw.pop() to retrieve input arguments from kw.
        # --- This allows testing for any unused arguments and raising an error if found.
        pass
