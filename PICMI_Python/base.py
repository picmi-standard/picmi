"""base code for the PICMI standard
"""

codename = None

# --- The list of supported codes is needed to allow checking for bad arguments.
supported_codes = ['warp', 'warpx', 'fbpic']

def register_codename(_codename):
    global codename
    codename = _codename

class _ClassWithInit(object):
    def handle_init(self, kw):
        # --- Grab all keywords for the current code.
        # --- Arguments for other supported codes are ignored.
        # --- If there is anything left over, it is an error.
        codekw = {}
        for k,v in kw.copy().items():
            code = k.split('_')[0]
            if code == codename:
                codekw[k] = v
                kw.pop(k)
            elif code in supported_codes:
                kw.pop(k)

        if kw:
            raise TypeError('Unexpected keyword argument: "%s"'%list(kw))

        # --- It is expected that init strips accepted keywords from codekw.
        self.init(codekw)

        if codekw:
            raise TypeError("Unexpected keyword argument for %s: '%s'"%(codename, list(codekw)))

    def init(self, kw):
        # --- The implementation of this routine should use kw.pop() to retrieve input arguments from kw.
        # --- This allows testing for any unused arguments and raising an error if found.
        pass
