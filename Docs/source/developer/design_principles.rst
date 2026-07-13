Design Principles
=================

This document describes fundamental design choices of PICMI.

Main Structure
--------------

PICMI defines a Python interface to specify a PIC input parameter set.
It is not a pure schema, i. e. can do computation.

In practical terms, PICMI defines a set of classes which implementations
must inherit from with well-defined names. These classes’ main purpose
is to provide a common user experience.

Semantics
---------

PICMI is (aims to be) well-defined enough to compare codes.

Scope
-----

PICMI aims to define all parameters for most PIC scenarios. Obscure
parameters may be implemented by codes when needed.

   Rule of thumb: If multiple codes implement a feature, put it into
   PICMI.

Grouping Principles
-------------------

Parameters defining physics are separated from parameters defining
numerics (i. e. how these physics are represented in the simulation).

Parameters should be optional if possible.

Parameters are grouped conveniently for the user: Classes may contain more
than the *minimal required set of parameters*, or in other terms:
Forming sub-groups of parameters is only done when **convenient**, not
mereley because sub-grouping would be possible.

One PICMI object *should* represent one physical thing.

Computations inside PICMI
-------------------------

PICMI defines a parameter structure. Convenient computations are
performed automatically.

   "Computations" here are filling in different representations of the
   same input data.

Redundancies are allowed, though should be limited to different representations of a physical concept (e.g. laser a0/e0).
Inconsistencies must be caught by checks when invoked explicitly.

Python Features
---------------

PICMI is implemented using Python. This section addresses some practical
approaches to the behavior and treatment of PICMI objects.

PICMI Object Lifecycle
~~~~~~~~~~~~~~~~~~~~~~

0. Parameters are collected/loaded/created by the user.
1. PICMI object is constructed, **all** parameters are passed.
   Automated computations are invoked.
   Their results and the original parameters are stored to be retrievable by the code.
   Note that the variables are exposed to the implementing code, not necessarily to the user.
2. The implementing code extracts these parameters.

PICMI objects mainly hold data and are not designed to be modified
after they have been initialized with their content.

Note that even if not recommended direct access to member variables
of objects is still possible.
Accessing the member variables directly (both reading and writing) results in **undefined behavior**.

Implementing codes may invoke checks on PICMI objects to ensure their integrity.

Types
~~~~~

For simplicity **strong typing** is used, i.e. variables are checked
against a finite list of type specifications. However, this finite list
of type specifications is kept **permissive**, including (1) general
types (e.g. iterable instead of "list"), and (2) common library types
(e.g. numpy types).

Identifying Objects
~~~~~~~~~~~~~~~~~~~

Two PICMI objects are considered to be equal (and hence should yield the
same simulation initialization) if they fullfill the Python equals
operator (``__eq__()``, invoked by ``==``).

PICMI reserves the option to redefine ``__eq__()`` to follow either of
these approaches:

- Two PICMI objects are considered equal if all their attributes are equal (dataclass approach).
- Two PICMI objects are considered equal if they point to the same object,
  i.e. the same memory region (python ``is`` operator, python vanilla for custom objects).

Friendly towards Implementations
--------------------------------

It should be *easy* to implement PICMI.

Extensible
~~~~~~~~~~

PICMI passes through all arguments (including codespecific variables)
which can be handled by implementations.

It is not a library, i. e. it does not provide functionality to be used
by the implementation.

Safety
~~~~~~

PICMI aims to provide a well defined interface to implementing codes.
PICMI guarantees the integrity of its structure by providing ``check()``
methods, which can be called explicitly.

PICMI is **predictable** in the sense that any variable provided by the
user is passed through to the implementation using **the same name**.

Implementations are not Bound by PICMI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementations may extend PICMI at arbitrary points. PICMI provides
simple interfaces for that. Implementations may only implement PICMI
partially.
