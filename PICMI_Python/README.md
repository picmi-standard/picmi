# PICMI: Python implementation

The *Particle-In-Cell Modeling Interface* (**PICMI**) standard establishes conventions for the naming and structuring of input files for Particle-In-Cell simulations.

The goal of the standard is to propose a set (or dictionary) of names and definitions that can be used in simulations, with as little changes as possible between codes.

This contains the base classes for the implementation of the PICMI standard. In order to install this package, run:
  ```
  python setup.py install
  ```
  The latest release of this package is also available via `pip`:
  ```
  pip install picmistandard
  ```

The Python classes provided here should be the bases classes of the code specific implementation. For each class, the implementation
should take the inputs and convert them to the appropriate inputs for the specific code.

**Warning:**
The standard is still evolving at this point. The most basic components of the standard have been defined and are being refined. New components are being added.

