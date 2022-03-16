General Class Design
====================

All PICMI classes inherit from ``_ClassWithInit``, which makes a class defined by purely attributes exhibit *the typical PICMI behavior*:

.. code-block:: python

   >>> class MyClass(_ClassWithInit):
   ...     mandatory_attr: typing.Any
   ...     name: str = ""
   ...     optional: typing.Optional[int] = 3
   ...     _protected = 1
   >>> my_object = MyClass(mandatory_attr=[],
   ...                     name="any string")
   >>> [] == my_object.mandatory_attr
   True
   >>> "any string" == my_object.name
   True
   >>> 3 == my_object.optional
   True
   >>> my_object.optional = None
   >>> my_object.check() # <- silently passes
   >>> my_object.optional = "not an integer"
   >>> my_object.check() # <- now throws
   [...TypeError stack trace omitted...]


General Functionality
---------------------
When a class has a set attributes, which may have a `PEP 484 type annotation <https://peps.python.org/pep-0484/>`_, and a default value, a constructor is provided which behaves as follows:

For all arguments (``kwargs``) check:

- **allow** if:
  - attribute is defined, or
  - attribute starts with the current code name, or
  - attribute starts with a known code name
- **reject** if:
  - attribute starts with a ``_`` (is private or protected), or
  - is unkown (and not one of the exceptions above)

If a type annotation is given, the value is checked against it, otherwise all types are allowed.

If a type has **no default value** it is considered **mandatory**:
The constructure will abort operation if it is not given.

Additionally, an implementing class may implement a method ``_check(self) -> None``, which can perform additional checks.

All checks (type checks and custom checks by ``_check()``) can be manually invoked by calling ``check()``.

After construction ``check()`` will automatically be called.

An implementation should call ``check()`` just before it begins it operation too.

Mandatory and Optional Attributes
---------------------------------

Unless a default value is given, attributes are considered mandatory.
Note that the default must conform to the given type specification -- to set a default to ``None``, the type specification must allow it.

Type Checking
-------------

Type checking is performed only for attributes, and delegated to the library `typeguard <https://typeguard.readthedocs.io/en/latest/>`_.
Methods must be checked by other means (e.g. by using typeguard ``@typechecked`` annotation).

Other Checks
------------

The method ``check()`` provides type checking for the constructor and other objects.

- ``check()``: **external** interface, performs typechecks and calls ``_check()``. **DO NOT OVERWRITE**
- ``_check()``: **custom** hook, optional, automatically called from wrapper ``check()``. **DO OVERWRITE**

To implement custom checks overwrite ``_check(self) -> None``, which will be called from the wrapper ``check()`` after all typechecks have passed.

``_check()`` must raise when it encounters an error, passing silently will consider the check to have succeeded.

Note that owned objects will not be checked by default, consider calling their respective ``check()`` inside the custom ``_check()``.

Full Reference
--------------

.. autoclass:: picmistandard.base._ClassWithInit
   :members: __init__, _check, check
