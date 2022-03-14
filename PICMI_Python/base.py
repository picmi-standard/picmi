"""base code for the PICMI standard
"""
import inspect
import warnings

codename = None

# --- The list of supported codes is needed to allow checking for bad arguments.
supported_codes = ['warp', 'warpx', 'fbpic']

def register_codename(_codename):
    """This must be called by the implementing code, passing in the code name"""
    global codename
    codename = _codename

# --- This needs to be set by the implementing package (by calling register_constants).
# --- It allows constants to be used within the picmi interface, with the constants
# --- defined in the implementation.
_implementation_constants = None

def register_constants(implementation_constants):
    """This must be called by the implementing code, passing in the constans object"""
    global _implementation_constants
    _implementation_constants = implementation_constants

def _get_constants():
    return _implementation_constants

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

    def _check_unsupported_argument(self, arg_name, message=None, raise_error=False):
        """Raise a warning or exception if an unsupported argument was specified by the user
        - arg_name: The name of the unsupported argument (string)
        - message: Information to include in the warning/error message (string)
        - raise_error: If False (the default), raise a warning. If true, raise an exception
                       (which interrupts the code).

        Implementation note: This should be called in the "init" method of the
        implementing class for each unsupported argument. For example, for the
        'density_scale' argument of Species:

            self._check_unsupported_argument(
                'density_scale',
                message='My code can not handle a density_scale')

        """

        # This compares the value of the parameter with the dault value in the __init__ method.
        # If they differ, this means that the user supplied a value, so a warning or error is raised.
        signature = inspect.signature(self.__init__)
        default_value = signature.parameters[arg_name].default
        if not (getattr(self, arg_name) == default_value):
            message = f'{self.__name__}: For argument {arg_name} is not supported.'
            if message is not None:
                message += f' {message}'
            if raise_error:
                raise Exception(message)
            else:
                warnings.warn(message)

    def _unsupported_value(self, arg_name, message='', raise_error=True):
        """Raise a warning or exception for argument with an unsupported value.
        - arg_name: The name of the argument with an unsupported value (string)
        - message: Information to include in the warning/error message (string)
        - raise_error: If False (the default), raise a warning. If true, raise an exception
                       (which interrupts the code).

        Implementation note: This should be called when the implementing code handles
        the input arguments. For example, for 'method' in Species:

            if self.method not in ['Boris', 'Li']:
                self._unsupported_value(
                    'method',
                    message='My code only supports Boris and Li')

        """
        message = f'{self.__name__}: For argument {arg_name}, the value {getattr(self, arg_name)} is not supported.'
        if message is not None:
            message += f' {message}'
        if raise_error:
            raise Exception(message)
        else:
            warnings.warn(message)

    def _check_deprecated_argument(self, arg_name, message=None, raise_error=False):
        """Raise a warning or exception if a deprecated argument was specified by the user
        - arg_name: The name of the deprecated argument (string)
        - message: Information to include in the warning/error message (string)
        - raise_error: If False (the default), raise a warning. If true, raise an exception
                       (which interrupts the code).

        Implementation note: This should be called within PICMI in the "__init__" method of classes
        for each deprecated argument. This assumes that the argument is still included in the
        argument list of __init__ as a transition until it is removed.
        For example, if the 'density_scale' argument of Species was to be deprecated:

            self._check_deprecated_argument(
                'density_scale',
                message='This argument is no longer needed')

        """

        # This compares the value of the parameter with the dault value in the __init__ method.
        # If they differ, this means that the user supplied a value, so a warning or error is raised.
        signature = inspect.signature(self.__init__)
        default_value = signature.parameters[arg_name].default
        if not (getattr(self, arg_name) == default_value):
            message = f'{self.__name__}: For argument {arg_name} is not supported.'
            if message is not None:
                message += f' {message}'
            if raise_error:
                raise Exception(message)
            else:
                warnings.warn(message)

