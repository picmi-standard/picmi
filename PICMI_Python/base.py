"""base code for the PICMI standard
"""

import typing
import typeguard

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


class _ClassWithInit:
    """
    Use args from constructor as attributes

    Non-given attributes are left untouched (i.e. at their default).

    Attributes can be marked as mandatory by adding a type annotation.
    Type annotations are enforced by check(), i.e. normal assignments for
    attributes with wrong types still work.

    Arguments that are prefixed with codename and underscore _ will be
    accepted, as well as equivalent prefixes for other supported codes.
    All other arguments will be rejected.
    """

    def __check_arg_valid(self, arg_name: str) -> None:
        """
        check if arg_name is acceptable for attr

        If ok silently pass, else raise.

        An arg is valid if:
        - has a default
        - has a type annotation (i.e. is mandatory)
        - is prefixed with codename and underscore _
        - is prefixed with any supported codename and underscore _
        - is equal to codename/any supported codename
        - does *NOT* begin with underscore _
        """
        assert codename is not None
        self_type = type(self)

        if arg_name.startswith("_"):
            raise NameError(
                f"protected/private attribute may NOT be accessed: {arg_name}")

        if arg_name in self_type.__dict__:
            # has default value i.e. is defined
            return

        if arg_name in typing.get_type_hints(self_type):
            # mandatory arg (has no default)
            return

        # check prefix:
        prefix = arg_name.split("_")[0]
        if prefix == codename:
            return

        if prefix in supported_codes:
            return

        # arg name is not in allowed sets -> raise
        raise NameError(f"unkown argument: {arg_name}")

    def __get_mandatory_attrs(self) -> typing.Set[str]:
        """
        Retrieve list of mandatory attrs

        Attributes are considered mandatory if they exist (which they
        syntactically only can if they have a *type annotation*), but no
        default value.
        """
        self_type = type(self)
        has_type_annotion = typing.get_type_hints(type(self)).keys()

        # ignore those with default values
        return set(has_type_annotion - self_type.__dict__.keys())

    def check(self) -> None:
        """
        checks self, raises on error, passes silently if okay

        ! MUST NOT BE OVERWRITTEN !
        -> overwrite _check() (note the leading _) instead

        Performs the following checks:
        1. The type of all attributes are checked against their type
           annotations (if present). On error a TypeError is raised.
        2. If (and only if) the typechecks passed the custom-overwriteable
           method _check() is called.

        Will be called inside of __init__(), and should be called before any
        work on the data is performed.
        """
        for arg_name, annotation in typing.get_type_hints(type(self)).items():
            # perform type check: will raise on error
            typeguard.check_type(arg_name,
                                 getattr(self, arg_name),
                                 annotation)

        self._check()

    def __check_types_of_defaults(self) -> None:
        """
        run typechecks for defaults

        Ensures that defaults for attributes conform to their type annotations.
        Passes silently if okay, else raises type error.
        """
        for arg_name, annotation in typing.get_type_hints(type(self)).items():
            if arg_name not in type(self).__dict__:
                # there is no default -> skip loop iteration
                continue
            default_value = type(self).__dict__[arg_name]
            # perform type check: will raise on error
            try:
                typeguard.check_type(arg_name, default_value, annotation)
            except TypeError:
                # replace by custom message
                actual_type = type(default_value)
                raise TypeError(f"default value for {arg_name} must be of "
                                "type {annotation}, but got {actual_type}")

    def _check(self) -> None:
        """
        run checks that are not typechecks

        Should be overwritten by child class.

        Will be called from check(), and thereby from __init__().

        When this method passes it guarantees that self is conforming to PICMI.
        This includes all children (should there be any) -- their check()
        method should be called from this method.

        When it is not overwritten, it is replaced by this empty parent method.
        """
        # parent implementation: just pass
        pass

    def __init__(self, **kw):
        """
        parse kw and set class attributes accordingly

        See class docstring for detailed description.

        ! MUST NOT BE OVERWRITTEN !
        (Constructors MUST NOT exhibit unpredictable behavior==behavior
        different from the one specified here.)
        """
        mandatory_missing = self.__get_mandatory_attrs() - kw.keys()
        if 0 != len(mandatory_missing):
            raise RuntimeError(
                "mandatory attributes are missing: {}"
                .format(", ".join(mandatory_missing)))

        self.__check_types_of_defaults()

        for name, value in kw.items():
            self.__check_arg_valid(name)
            setattr(self, name, value)

        # perform self-check -> will alert for invalid params
        self.check()
