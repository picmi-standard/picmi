"""base code for the PICMI standard
"""
from __future__ import annotations
from typing import Any
import sys

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

import warnings

from pydantic import BaseModel, ConfigDict, Field, PrivateAttr, model_validator
from pydantic._internal._model_construction import ModelMetaclass


codename: str | None = None

# --- The list of supported codes is needed to allow checking for bad arguments.
supported_codes = ['warp', 'warpx', 'fbpic']

def register_codename(_codename: str) -> None:
    """This must be called by the implementing code, passing in the code name"""
    global codename
    codename = _codename

# --- This needs to be set by the implementing package (by calling register_constants).
# --- It allows constants to be used within the picmi interface, with the constants
# --- defined in the implementation.
_implementation_constants: Any = None

def register_constants(implementation_constants: Any) -> None:
    """This must be called by the implementing code, passing in the constans object

    Parameters
    ----------
    implementation_constants: python object
        The object must have as attributes the physical constants
    """
    global _implementation_constants
    _implementation_constants = implementation_constants

def _get_constants() -> Any:
    return _implementation_constants


class _DocumentedMetaClass(ModelMetaclass):
    """This is used as a metaclass that combines the __doc__ of the picmistandard base and of the implementation"""
    def __new__(mcs, name, bases, namespace, **kwargs):
        # For _ClassWithInit itself, override BaseModel's docstring to remove "Usage docs" link
        if name == '_ClassWithInit' and 'BaseModel' in [b.__name__ for b in bases if hasattr(b, '__name__')]:
            # Use only our custom docstring, not BaseModel's
            namespace['__doc__'] = namespace.get('__doc__', '')
        
        # "if bases" skips this for the _ClassWithInit (which has no bases)
        # "if bases[0].__doc__ is not None" skips this for the picmistandard classes since their bases[0] (i.e. _ClassWithInit)
        # has no __doc__.
        elif bases and bases[0].__doc__ is not None:
            implementation_doc = namespace.get('__doc__', '')
            if implementation_doc:
                # The format of the added string is intentional.
                # The double return "\n\n" is needed to start a new section in the documentation.
                # Then the four spaces matches the standard level of indentation for doc strings
                # (assuming PEP8 formatting).
                # The final return "\n" assumes that the implementation doc string begins with a return,
                # i.e. a line with only three quotes, """.
                namespace['__doc__'] = bases[0].__doc__ + """\n\n    Implementation specific documentation\n""" + implementation_doc
            else:
                namespace['__doc__'] = bases[0].__doc__
        return super().__new__(mcs, name, bases, namespace, **kwargs)


class _ClassWithInit(BaseModel, metaclass=_DocumentedMetaClass):
    """
    Base class for all PICMI classes using Pydantic for validation and extensibility.
    
    This class allows code-specific extensions (e.g., warpx_* kwargs) via Pydantic's
    extra fields mechanism while maintaining type safety for standard attributes.
    """
    
    model_config = ConfigDict(
        # Allow extra fields for code-specific extensions
        extra='allow',
        # Validate assignment to catch errors early
        validate_assignment=True,
        # Use enum values instead of enum names
        use_enum_values=True,
        # Allow arbitrary types for flexibility (e.g., grid objects, distribution objects)
        arbitrary_types_allowed=True,
    )
    
    # Store code-specific kwargs separately for backward compatibility
    _code_kwargs: dict[str, Any] = PrivateAttr(default_factory=dict)
    
    @model_validator(mode='before')
    @classmethod
    def _handle_code_specific_kwargs(cls, data: Any) -> Any:
        """
        Process code-specific kwargs (e.g., warpx_*) before validation.
        This maintains backward compatibility with the old handle_init pattern.
        """
        if not isinstance(data, dict):
            return data
            
        # Make a copy to avoid mutating the input
        data = dict(data)
        
        # Separate code-specific kwargs from standard kwargs
        codekw: dict[str, Any] = {}
        standard_kwargs: dict[str, Any] = {}
        
        for k, v in data.items():
            if '_' in k:
                code = k.split('_')[0]
                if code == codename:
                    codekw[k] = v
                elif code in supported_codes:
                    # Ignore kwargs for other supported codes
                    pass
                else:
                    standard_kwargs[k] = v
            else:
                standard_kwargs[k] = v
        
        # Store code-specific kwargs for later processing
        if codekw:
            standard_kwargs['_code_kwargs'] = codekw
        
        return standard_kwargs
    
    @model_validator(mode='after')
    def _process_code_kwargs(self) -> 'Self':
        """
        Process code-specific kwargs after validation.
        This can be overridden in subclasses to handle code-specific attributes.
        """
        if hasattr(self, '_code_kwargs') and self._code_kwargs:
            # Call the init method for backward compatibility
            self.init(self._code_kwargs)
            
            # Check for unused code-specific kwargs
            if self._code_kwargs:
                raise TypeError(
                    f"Unexpected keyword argument for {codename}: '{list(self._code_kwargs)}'"
                )
        return self
    
    def init(self, kw: dict[str, Any]) -> None:
        """
        The implementation of this routine should use kw.pop() to retrieve input 
        arguments from kw. This allows testing for any unused arguments and raising 
        an error if found.
        
        Subclasses can override this method to handle code-specific arguments.
        """
        pass
    
    def handle_init(self, kw: dict[str, Any]) -> None:
        """
        Legacy method for backward compatibility.
        In Pydantic-based classes, this is handled automatically during model validation.
        """
        # This method is kept for backward compatibility but is largely unused
        # in the Pydantic implementation as validation happens automatically
        pass
    
    def _check_unsupported_argument(self, arg_name: str, message: str | None = None, raise_error: bool = False) -> None:
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

        Implementation note: This should be called in the "init" method of the
        implementing class for each unsupported argument. For example, for the
        'density_scale' argument of Species:

            self._check_unsupported_argument(
                'density_scale',
                message='My code can not handle a density_scale')

        """
        # Get the field info from the Pydantic model
        if arg_name in self.model_fields:
            field_info = self.model_fields[arg_name]
            default_value = field_info.default
            if default_value is not ...:
                current_value = getattr(self, arg_name, None)
                if not (current_value == default_value):
                    full_message = f'{self.__class__.__name__}: For argument {arg_name} is not supported.'
                    if message is not None:
                        full_message += f' {message}'
                    if raise_error:
                        raise Exception(full_message)
                    else:
                        warnings.warn(full_message)
    
    def _unsupported_value(self, arg_name: str, message: str = '', raise_error: bool = True) -> None:
        """Raise a warning or exception for argument with an unsupported value.

        Parameters
        ----------
        arg_name: string
            The name of the argument with an unsupported value

        message: string
            Information to include in the warning/error message

        raise_error: bool
            If False (the default), raise a warning. If true, raise an exception
            (which interrupts the code).

        Implementation note: This should be called when the implementing code handles
        the input arguments. For example, for 'method' in Species:

            if self.method not in ['Boris', 'Li']:
                self._unsupported_value(
                    'method',
                    message='My code only supports Boris and Li')

        """
        current_value = getattr(self, arg_name, None)
        full_message = f'{self.__class__.__name__}: For argument {arg_name}, the value {current_value} is not supported.'
        if message:
            full_message += f' {message}'
        if raise_error:
            raise Exception(full_message)
        else:
            warnings.warn(full_message)
    
    def _check_deprecated_argument(self, arg_name: str, message: str | None = None, raise_error: bool = False) -> None:
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

        Implementation note: This should be called within PICMI in the "__init__" method of classes
        for each deprecated argument. This assumes that the argument is still included in the
        argument list of __init__ as a transition until it is removed.
        For example, if the 'density_scale' argument of Species was to be deprecated:

            self._check_deprecated_argument(
                'density_scale',
                message='This argument is no longer needed')

        """
        if arg_name in self.model_fields:
            field_info = self.model_fields[arg_name]
            default_value = field_info.default
            if default_value is not ...:
                current_value = getattr(self, arg_name, None)
                if not (current_value == default_value):
                    full_message = f'{self.__class__.__name__}: For argument {arg_name} is deprecated.'
                    if message is not None:
                        full_message += f' {message}'
                    if raise_error:
                        raise Exception(full_message)
                    else:
                        warnings.warn(full_message)

