Implementing the PICMI standard in an existing code
===================================================

.. warning::

   This section is currently in development.

In order to implement the PICMI standard, a given code should go through
the following steps:

- Define a **Python package** (ideally, whose name is the name of the code itself), which contains a module ``picmi``:

    The ``picmi`` module should be importable with the following syntax:

    ::

        from <python_package> import picmi

    where ``<python_package>`` should be replaced by the name of the code.

- In the ``picmi`` module, define the variable ``codename`` as a string containing the name of the code,
  and register it with ``picmistandard.register_codename``.
  Keyword arguments that start with the name of another supported code (``picmistandard.base.supported_codes``) are then ignored,
  so that a script can carry the arguments of several codes.

- Define a class ``constants`` with the constants described in :doc:`../standard/constants`,
  and register it with ``picmistandard.register_constants``.
  The standard classes use it, e.g., to compute the field amplitude of a laser from its normalized vector potential.

- Define the **PICMI classes** that a user of this code would typically use, among those defined in :doc:`../standard/standard`.

    This is done by deriving a class from the corresponding class in the package ``picmistandard``.
    The PICMI classes are `pydantic <https://docs.pydantic.dev>`__ models:

    - Code-specific parameters are fields of the derived class.
      Users give them with the name of the code as a prefix (``<codename>_``), which is the alias of the field.
      They can read and assign them under that name, too, e.g., ``simulation.mycode_load_balance_interval``,
      while the code itself uses the name of the field, e.g., ``self.load_balance_interval``.
    - A different default of a standard parameter is given by declaring the field again with this default.
    - Code-specific initialization goes into ``model_post_init``, which is called after the parameters are validated.
    - The state of an object that is not a parameter is kept in private attributes, whose names start with ``_``.

    For instance:

    ::

        import picmistandard
        from pydantic import Field, PrivateAttr

        codename = "mycode"
        picmistandard.register_codename(codename)


        class constants:
            c = 299792458.0
            ep0 = 8.8541878188e-12
            mu0 = 1.2566370612685e-06
            q_e = 1.602176634e-19
            m_e = 9.1093837139e-31
            m_p = 1.67262192595e-27


        picmistandard.register_constants(constants)


        class Simulation(picmistandard.PICMI_Simulation):
            # a code-specific parameter, given as ``mycode_load_balance_interval``
            load_balance_interval: int | None = Field(
                default=None,
                alias="mycode_load_balance_interval",
                description="Number of steps between load balancing",
            )

            # code-specific state, which is not a parameter
            _initialized_space_charge_sources: bool = PrivateAttr(default=False)

            def model_post_init(self, context):
                super().model_post_init(context)
                # code-specific initialization
                ...

- Set ``picmistandard.PICMI_MultiSpecies.Species_class`` to the species class of the code,
  which ``MultiSpecies`` creates its species with.

- Derive classes that have no counterpart in the standard, e.g., additional field solvers or diagnostics,
  from the extension class of their kind (see :doc:`../standard/extensions`).

.. note::

    For concrete examples on how to implement the PICMI standard, see the implementation in:

        - `WarpX <https://github.com/BLAST-WarpX/warpx/blob/development/Python/pywarpx/picmi.py>`__
        - `Warp <https://bitbucket.org/berkeleylab/warp/src/master/scripts/picmi.py>`__
        - `FBPIC <https://github.com/fbpic/fbpic/tree/dev/fbpic/picmi>`__
        - `PIConGPU <https://github.com/ComputationalRadiationPhysics/picongpu/tree/dev/lib/python/picongpu/picmi>`__
