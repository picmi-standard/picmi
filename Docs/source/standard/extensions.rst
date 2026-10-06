Code-specific extensions
========================

Implementing codes can provide classes that have no counterpart in the standard, e.g., additional field solvers or diagnostics.
They derive these classes from the base class of the respective kind, which the classes of the standard derive from as well, so that the objects are accepted by the PICMI classes where objects of that kind are expected.
Each of these base classes is part of the type alias of its kind (see :doc:`types`), e.g., :py:class:`~picmistandard.PICMI_Solver` in :py:data:`~picmistandard.PICMI_AnySolver`.
Classes that are only used by other code-specific classes derive from ``PICMI_Extension`` directly.

.. autopydantic_model:: picmistandard.PICMI_Extension
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Grid
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Solver
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Distribution
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Layout
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Laser
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_LaserInjection
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_AppliedField
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Diagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_Interaction
    :inherited-members: BaseModel

Classes with analytic expressions, whose parameters are given as additional keyword arguments, derive from ``PICMI_ExpressionParameters``, too.

.. autopydantic_model:: picmistandard.PICMI_ExpressionParameters
    :inherited-members: BaseModel
