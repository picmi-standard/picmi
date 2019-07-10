Implementing the PICMI standard in an existing code
===================================================

.. warning::

   This section is currently in development.

In order to implement the PICMI standard, a given code should go through
the following step:

- Define a **Python package** (ideally, whose name is the name of the code itself), which contains a module ``picmi``:

    The ``picmi`` module should be importable with the following syntax:

    ::

        from <python_package> import picmi

    where ``<python_package>`` should be replaced by the name of the code.

- Define the variable ``picmi.codename`` as a string containing the name of the code.

- Create a sub-module ``picmi.constants``, that defines the constants described in :doc:`../standard/constants.rst`.

- Define the **PICMI classes** that a user of this code would typically us, among those defined in :doc:`../standard`.

    This is should be done by defining a subclass of the corresponding class in the package ``picmistandard``,
    and defining its ``init`` method. For instance:

    ::

        import picmistandard

        class Simulation( picmistandard.PICMI_Simulation ):

            def init(self, kw):

                # Code-specific initialization
                self.initialized_space_charge_sources = False
