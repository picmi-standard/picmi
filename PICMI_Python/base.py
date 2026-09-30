"""base code for the PICMI standard
"""
import contextlib
import functools
import json
import numbers
import re
import threading
from itertools import repeat
import warnings
from typing import Annotated, ClassVar, Self
from pydantic import BeforeValidator, model_serializer, model_validator, BaseModel, ConfigDict
from typing_extensions import TypeAliasType


# Tracks which (instance, validator) pairs are currently executing, per thread, so that a
# mode="after" validator that assigns to ``self`` does not re-enter itself when
# ``validate_assignment=True`` re-validates each of those assignments.
_validators_in_progress = threading.local()


def resolve_once(validator):
    """Make a ``mode="after"`` model validator safe under ``validate_assignment=True``.

    With assignment validation enabled, every ``self.x = ...`` performed inside an
    after-validator triggers a re-validation, which re-runs the very same validator and
    would recurse without bound. This wrapper turns a re-entrant call *on the same
    instance* into a no-op (it returns ``self`` unchanged), so the validator's own
    derived-field assignments do not re-execute its body. The outermost call still runs in
    full, so derived fields are computed/resolved exactly once per validation.

    The fields that the validator derives are not marked as set (``model_fields_set``), so
    that dumps with ``exclude_unset=True`` contain only the parameters that were given.

    Apply it *under* ``@model_validator(mode="after")``::

        @model_validator(mode="after")
        @resolve_once
        def _resolve(self) -> Self:
            ...
    """
    @functools.wraps(validator)
    def wrapper(self):
        active = _validators_in_progress.__dict__.setdefault("markers", set())
        marker = (id(self), validator)
        if marker in active:
            return self
        active.add(marker)
        fields_set = set(self.__pydantic_fields_set__)
        try:
            return validator(self)
        finally:
            active.discard(marker)
            object.__setattr__(self, "__pydantic_fields_set__", fields_set)
    return wrapper

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
    """This must be called by the implementing code, passing in the constans object

    Parameters
    ----------
    implementation_constants: python object
        The object must have as attributes the physical constants
    """
    global _implementation_constants
    _implementation_constants = implementation_constants

def _get_constants():
    return _implementation_constants


class _DocumentedModelMetaClass(type(BaseModel)):
    """Metaclass that combines the __doc__ of the picmistandard base and of the implementation.

    Downstream codes (e.g. WarpX) can extend the documentation of a PICMI class simply by
    adding a docstring to their subclass. It derives from pydantic's metaclass
    (``type(BaseModel)`` is ``ModelMetaclass``) so that it composes with ``BaseModel``.
    """
    def __new__(mcs, name, bases, namespace, **kwargs):
        # Skip the infrastructure base itself (its only base is BaseModel), any class whose
        # first base carries no docstring (e.g. _PICMIModel), and the base classes whose
        # docstrings describe a mechanism (e.g. the extensions) rather than the derived class.
        if (
            bases
            and bases[0] is not BaseModel
            and bases[0].__doc__ is not None
            and not bases[0].__dict__.get("__picmi_doc_not_inherited__", False)
        ):
            implementation_doc = namespace.get('__doc__', '')
            if implementation_doc:
                # The double return "\n\n" separates the picmistandard docstring from the
                # implementation-specific one, starting a new paragraph in the documentation.
                namespace['__doc__'] = bases[0].__doc__ + "\n\n" + implementation_doc
            else:
                namespace['__doc__'] = bases[0].__doc__
        return super().__new__(mcs, name, bases, namespace, **kwargs)


# --- Serialized PICMI objects carry the class they were dumped from under this key, so
# --- that loading them restores the same (e.g. code-specific) class, also when nested in
# --- another PICMI object or in a field that is not typed with a specific class.
PICMI_CLASS_KEY = "picmi_class"

# --- All pydantic-based PICMI classes (of the standard and of the implementing codes),
# --- by their _picmi_class_name. Only these classes are instantiated when loading data.
_picmi_classes = {}


def _picmi_class_name(cls):
    return f"{cls.__module__}.{cls.__qualname__}"


def _registered_picmi_class(name):
    try:
        return _picmi_classes[name]
    except KeyError:
        raise ValueError(
            f"Unknown {PICMI_CLASS_KEY} '{name}'. Import the module that defines this class before loading the data."
        ) from None


def load(data):
    """Load a PICMI object as the class that it was dumped from

    Parameters
    ----------
    data: str, bytes or dict
        The JSON of the object (from ``model_dump_json``) or its dictionary (from ``model_dump``).

    The module that defines the class, e.g., of the implementing code, must be imported before.
    """
    parsed = json.loads(data) if isinstance(data, (str, bytes, bytearray)) else data
    if not isinstance(parsed, dict) or PICMI_CLASS_KEY not in parsed:
        raise ValueError(f"The data does not record the class of a PICMI object ({PICMI_CLASS_KEY}).")
    picmi_class = _registered_picmi_class(parsed[PICMI_CLASS_KEY])
    if parsed is data:
        return picmi_class.model_validate(data)
    return picmi_class.model_validate_json(data)


def _instantiate_picmi_objects(value):
    """Turn (possibly nested in lists or tuples) dictionaries written by ``model_dump`` into
    instances of the PICMI class recorded in them."""
    if isinstance(value, dict) and PICMI_CLASS_KEY in value:
        data = dict(value)
        # Remove the class marker already here: pydantic passes the data of nested objects to
        # a custom __init__ of their class before the model validators run.
        picmi_class = _registered_picmi_class(data.pop(PICMI_CLASS_KEY))
        return picmi_class.model_validate(data)
    if isinstance(value, list):
        return [_instantiate_picmi_objects(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_instantiate_picmi_objects(item) for item in value)
    return value


class _PICMIModel(BaseModel, metaclass=_DocumentedModelMetaClass):
    # Shared configuration for all pydantic-based PICMI classes.
    # - ``extra="forbid"`` restores the old behaviour of raising on unexpected keyword
    #   arguments (pydantic's default silently ignores them).
    # - ``populate_by_name`` lets downstream codes expose extension inputs under a
    #   ``<code>_`` alias while keeping their internal attribute name.
    # - ``polymorphic_serialization`` serializes an object with the fields of its actual
    #   class, e.g., a downstream grid passed to a solver keeps its code-specific fields
    #   (by default, pydantic uses the fields of the annotated standard class).
    model_config = ConfigDict(
        populate_by_name=True,
        extra="forbid",
        validate_assignment=True,
        polymorphic_serialization=True,
    )

    # The field of each alias, e.g., the field "break_signals" of "warpx_break_signals",
    # so that the parameters of implementing codes can also be read and assigned under the
    # name that they are given with (pydantic aliases are keyword arguments only).
    _field_of_alias: ClassVar[dict[str, str]] = {}

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs):
        super().__pydantic_init_subclass__(**kwargs)
        _picmi_classes[_picmi_class_name(cls)] = cls
        fields = cls.model_fields
        cls._field_of_alias = {
            field.alias: name
            for name, field in fields.items()
            # a name that is a field itself always refers to that field
            if field.alias and field.alias != name and field.alias not in fields
        }

    def __getattr__(self, name):
        field = type(self)._field_of_alias.get(name)
        if field is not None:
            return getattr(self, field)
        return super().__getattr__(name)

    def __dir__(self):
        # so that the aliases are offered by interactive completion, too
        return [*super().__dir__(), *type(self)._field_of_alias]

    def __setattr__(self, name, value):
        # Pydantic applies an assignment before running the model validators, and keeps it
        # if they reject it. Restore the previous state in that case, so that a failed
        # assignment (including the assignments that validators make to derived fields)
        # leaves the object unchanged and valid.
        name = type(self)._field_of_alias.get(name, name)
        if name not in type(self).model_fields:
            return super().__setattr__(name, value)
        with self._atomic_update():
            super().__setattr__(name, value)

    @contextlib.contextmanager
    def _atomic_update(self):
        """Context in which several assignments either all succeed or leave the object unchanged"""
        previous_fields = dict(self.__dict__)
        previous_fields_set = set(self.__pydantic_fields_set__)
        previous_private = None if self.__pydantic_private__ is None else dict(self.__pydantic_private__)
        try:
            yield
        except Exception:
            object.__setattr__(self, "__dict__", previous_fields)
            object.__setattr__(self, "__pydantic_fields_set__", previous_fields_set)
            object.__setattr__(self, "__pydantic_private__", previous_private)
            raise

    @model_serializer(mode="wrap")
    def _serialize_with_picmi_class(self, handler):
        data = handler(self)
        if isinstance(data, dict):
            data = {PICMI_CLASS_KEY: _picmi_class_name(type(self)), **data}
        return data

    @model_validator(mode="before")
    @classmethod
    def _load_picmi_classes(cls, data):
        # Counterpart of _serialize_with_picmi_class: nested serialized PICMI objects are
        # loaded as the class they were dumped from, instead of the (standard) class a field
        # is annotated with, or a plain dictionary for fields that are not typed.
        if not isinstance(data, dict):
            return data
        data = dict(data)
        if PICMI_CLASS_KEY in data:
            dumped_class = _registered_picmi_class(data.pop(PICMI_CLASS_KEY))
            if not issubclass(cls, dumped_class):
                raise ValueError(
                    f"The data was dumped from {_picmi_class_name(dumped_class)} and cannot be loaded as {_picmi_class_name(cls)}."
                )
        return {key: _instantiate_picmi_objects(value) for key, value in data.items()}

    # PICMI objects are mutable handles to distinct entities of a simulation: two species
    # with identical parameters are still two species. Keep the identity-based equality and
    # hashing of the pre-pydantic classes (pydantic's default compares by value, which makes
    # mutable models unhashable), so that instances can be used as dictionary keys, e.g., to
    # give per-species diagnostic options. Compare ``model_dump()`` results to compare values.
    __eq__ = object.__eq__
    __hash__ = object.__hash__

    @model_validator(mode="before")
    @classmethod
    def _ignore_other_codes_arguments(cls, data):
        # Mirror the old handle_init() behaviour: keyword arguments prefixed with the name
        # of *another* supported code are silently ignored, so that a single PICMI input
        # script can carry code-specific arguments for several codes at once. Arguments
        # prefixed with the active codename (or otherwise unknown arguments) are left in
        # place and validated normally, so that genuine typos are still reported thanks to
        # ``extra="forbid"``.
        if not isinstance(data, dict):
            return data
        return {
            k: v for k, v in data.items()
            if not ((prefix := k.split('_')[0]) in supported_codes and prefix != codename)
        }

    def _is_given(self, arg_name):
        """Whether the argument has a value other than its default"""
        field = type(self).model_fields[arg_name]
        return getattr(self, arg_name) != field.get_default(call_default_factory=True)

    def _check_unsupported_argument(self, arg_name, message=None, raise_error=False):
        """Raise a warning or exception if an unsupported argument was specified by the user

        Parameters
        ----------
        arg_name: string
            The name of the unsupported argument

        message: string
            Information to include in the warning/error message

        raise_error: bool
            If False (the default), raise a warning. If true, raise an exception
            (which interrupts the code).

        Implementation note: This should be called by the implementing class, e.g., in its
        ``model_post_init``, for each unsupported argument. For example, for the
        'density_scale' argument of Species:

            self._check_unsupported_argument(
                'density_scale',
                message='My code can not handle a density_scale')

        """
        # A value that differs from the default means that the user supplied a value.
        if self._is_given(arg_name):
            full_message = f'{type(self).__name__}: The argument {arg_name} is not supported.'
            if message is not None:
                full_message += f' {message}'
            if raise_error:
                raise Exception(full_message)
            else:
                warnings.warn(full_message)

    def _unsupported_value(self, arg_name, message=None, raise_error=True):
        """Raise a warning or exception for argument with an unsupported value.

        Parameters
        ----------
        arg_name: string
            The name of the argument with an unsupported value

        message: string
            Information to include in the warning/error message

        raise_error: bool
            If True (the default), raise an exception (which interrupts the code).
            If False, raise a warning.

        Implementation note: This should be called when the implementing code handles
        the input arguments. For example, for 'method' in Species:

            if self.method not in ['Boris', 'Li']:
                self._unsupported_value(
                    'method',
                    message='My code only supports Boris and Li')

        """
        full_message = f'{type(self).__name__}: For argument {arg_name}, the value {getattr(self, arg_name)} is not supported.'
        if message is not None:
            full_message += f' {message}'
        if raise_error:
            raise Exception(full_message)
        else:
            warnings.warn(full_message)

    def _check_deprecated_argument(self, arg_name, message=None, raise_error=False):
        """Raise a warning or exception if a deprecated argument was specified by the user

        Parameters
        ----------
        arg_name: string
            The name of the deprecated argument

        message: string
            Information to include in the warning/error message

        raise_error: bool
            If False (the default), raise a warning. If true, raise an exception
            (which interrupts the code).

        Implementation note: This should be called within PICMI, e.g., in a model validator
        of the class, for each deprecated argument. This assumes that the argument is still a
        field of the class as a transition until it is removed.
        For example, if the 'density_scale' argument of Species was to be deprecated:

            self._check_deprecated_argument(
                'density_scale',
                message='This argument is no longer needed')

        """
        # A value that differs from the default means that the user supplied a value.
        if self._is_given(arg_name):
            full_message = f'{type(self).__name__}: The argument {arg_name} is deprecated.'
            if message is not None:
                full_message += f' {message}'
            if raise_error:
                raise Exception(full_message)
            else:
                warnings.warn(full_message)


def _as_expression(value):
    """Expressions are stored as strings without line breaks; numbers are accepted, too."""
    # numbers.Real includes NumPy scalars, e.g., numpy.float32, which are no Python floats
    if isinstance(value, numbers.Real) and not isinstance(value, bool):
        value = f'{value}'
    if isinstance(value, str):
        value = value.replace('\n', '')
    return value


# A named alias, so that the documentation shows "Expression" instead of the annotation.
Expression = TypeAliasType("Expression", Annotated[str, BeforeValidator(_as_expression)])
"""An analytic expression, given as a string (or a number)"""


def _expression_strings(value):
    """All strings in a value that is possibly nested in lists, tuples or dictionaries"""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from _expression_strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _expression_strings(item)


class PICMI_ExpressionParameters(_PICMIModel):
    """
    Base class of classes with analytic expressions, whose parameters are given as keyword arguments.

    Parameters used in the expressions can be given as additional keyword arguments. Those
    referenced in an expression are collected into the ``user_defined_kw`` field, which each
    derived class declares. Other unknown keyword arguments are still rejected. It is up to
    the implementing code to make sure that all parameters used in the expressions are
    defined.

    Derived classes list the names of their fields that hold expressions (strings, possibly
    nested in lists or dictionaries) in ``_expression_fields``. Keyword arguments named like a
    field are not collected, unless a derived class excludes that field in ``_parameter_names``.
    """

    __picmi_doc_not_inherited__ = True

    _expression_fields: ClassVar[tuple[str, ...]] = ()

    @classmethod
    def _parameter_names(cls, data):
        """The names (and aliases) of the fields that the given input data sets as parameters

        Keyword arguments with these names are never collected into ``user_defined_kw``. These
        are all fields, unless a derived class excludes fields that the input does not use, e.g.,
        fields that only apply to other types of an object.
        """
        fields = cls.model_fields
        return set(fields) | {field.alias for field in fields.values() if field.alias}

    @model_validator(mode="before")
    @classmethod
    def _collect_expression_parameters(cls, data):
        if not isinstance(data, dict):
            return data
        data = dict(data)

        fields = cls.model_fields
        known = cls._parameter_names(data)
        expressions = [
            expression
            for name in cls._expression_fields
            for key in {name, fields[name].alias or name}
            if key in data
            for expression in _expression_strings(data[key])
        ]

        user_defined_kw_key = fields["user_defined_kw"].alias or "user_defined_kw"
        if user_defined_kw_key not in data:
            user_defined_kw_key = "user_defined_kw"
        user_defined_kw = dict(data.get(user_defined_kw_key) or {})
        for key in list(data):
            if key in known:
                continue
            if any(re.search(r'\b%s\b' % re.escape(key), expression) for expression in expressions):
                user_defined_kw[key] = data.pop(key)
        # only set it if given or collected, so that the field is not marked as set otherwise
        if user_defined_kw or user_defined_kw_key in data:
            data[user_defined_kw_key] = user_defined_kw
        return data


# --------------------------------------------
# Base classes of code-specific extensions
# --------------------------------------------


class PICMI_Extension(_PICMIModel):
    """
    Base class of code-specific classes that have no counterpart in the PICMI standard.

    Implementing codes derive their own classes, e.g., additional field solvers or diagnostics,
    from one of the kind-specific extension classes below, so that these objects are accepted
    by the fields of the PICMI classes of that kind, e.g., a ``PICMI_SolverExtension`` as
    ``Simulation.solver``. Classes that are only used by fields of the implementing code itself
    derive from this class directly.

    Like all PICMI classes, extensions validate their parameters (unknown keyword arguments are
    rejected, assignments are validated) and are restored as their own class when loaded from
    a serialized simulation.

    Example:

    .. code-block:: python

        class HybridSolver(picmistandard.PICMI_SolverExtension):
            \"\"\"A code-specific field solver\"\"\"

            grid: picmistandard.PICMI_AnyGrid
            electron_temperature: float = Field(description="Electron temperature [eV]")

        simulation = Simulation(solver=HybridSolver(grid=grid, electron_temperature=10.0))
    """

    __picmi_doc_not_inherited__ = True


class PICMI_SolverExtension(PICMI_Extension):
    """
    Base class of code-specific field solvers, accepted as ``Simulation.solver``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_DistributionExtension(PICMI_Extension):
    """
    Base class of code-specific particle distributions, accepted as ``Species.initial_distribution``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_LayoutExtension(PICMI_Extension):
    """
    Base class of code-specific particle layouts, accepted as layout in ``Simulation.add_species``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_LaserExtension(PICMI_Extension):
    """
    Base class of code-specific laser profiles, accepted as laser in ``Simulation.add_laser``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_LaserInjectionExtension(PICMI_Extension):
    """
    Base class of code-specific laser injection methods, accepted as injection method in ``Simulation.add_laser``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_AppliedFieldExtension(PICMI_Extension):
    """
    Base class of code-specific applied fields, accepted by ``Simulation.add_applied_field``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_DiagnosticExtension(PICMI_Extension):
    """
    Base class of code-specific diagnostics, accepted by ``Simulation.add_diagnostic``.
    """

    __picmi_doc_not_inherited__ = True


class PICMI_InteractionExtension(PICMI_Extension):
    """
    Base class of code-specific interactions, accepted by ``Simulation.add_interaction`` and as ``Species.interactions``.
    """

    __picmi_doc_not_inherited__ = True


def broadcast_validation(values, condition, message="Condition not met."):
    if not all(condition(value) for value in values):
        raise ValueError(f"{message} You gave: {values}.")
    return values


def with_mutually_exclusive(*args, defaults=None, required=False):
    """Class decorator: at most one of the arguments differs from its default (None by default)

    With ``required=True``, exactly one of the arguments must be given.
    """
    def decorator(cls):
        def _mutually_exclusive(self) -> Self:
            # make sure we don't override previously implemented behaviour:
            parent_check = getattr(super(decorated, self), "_mutually_exclusive", None)
            if parent_check is not None:
                parent_check()
            if len(non_default := {arg: value for arg, default in zip(args, repeat(None) if defaults is None else defaults) if (value:=getattr(self, arg)) != default}) > 1:
                raise ValueError(f"The arguments {args} are mutually exclusive. You gave: {non_default=}.")
            if required and not non_default:
                raise ValueError(f"One of the arguments {args} must be given.")
            return self

        # Create the subclass under the name of the decorated class. A class statement would
        # name it after its local variable, which then leaks into the validation error titles,
        # the JSON schema and the documentation. The docstring is inherited through
        # _DocumentedModelMetaClass.
        decorated = type(cls)(
            cls.__name__,
            (cls,),
            {
                "__module__": cls.__module__,
                "__qualname__": cls.__qualname__,
                "_mutually_exclusive": model_validator(mode="after")(_mutually_exclusive),
            },
        )
        return decorated
    return decorator
