Diagnostics
===========

.. warning::

   This section is currently in development.

Standard PIC diagnostics
------------------------

.. autopydantic_model:: picmistandard.PICMI_ParticleDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_ParticleBoundaryScrapingDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_FieldDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_ElectrostaticFieldDiagnostic
    :inherited-members: BaseModel

Lab-frame diagnostics
---------------------

These diagnostics are used when running boosted-frame simulations.

.. autopydantic_model:: picmistandard.PICMI_LabFrameParticleDiagnostic
    :inherited-members: BaseModel

.. autopydantic_model:: picmistandard.PICMI_LabFrameFieldDiagnostic
    :inherited-members: BaseModel
