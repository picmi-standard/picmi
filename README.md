# PICMI

The *Particle-In-Cell Modeling Interface* (**PICMI**) standards establish conventions for the naming and structuring of input files for Particle-In-Cell and particle accelerator simulations.

The goal of the standard is to propose a set (or dictionary) of names and definitions that can be used in simulations, with as little changes as possible between codes.

Just as two persons do generally not use the exact same set of words from a dictionary (whatever the language), each code will use syntax (associated with unique functionalities of the code) that others will not. However, it is expected that, just like two persons will use and understand the same words and syntax for the most common topics, codes will share the same language for common definitions and tasks.

For example, it is expected that most PIC codes can share the same syntax to define a grid, a commonly-used field solver, a set of particles, etc. On the other hand, only a subset of the codes will offer, e.g. mesh refinement, field ionization, or code steering by the user. In this case, it is up to the implementation with each code to decide whether to ignore, raise a warning or an error, when an undefined statement is encountered. For more flexibility and robustness, conditional statements are available to the user to execute or read lines of the input script, based on the targeted code. 

The intent is for the standard to be agnostic of the language of implementation, which can then potentially be performed using JSON, Python, FORTRAN or other means.

**Warning:**
The standard is still evolving at this point. The most basic components of the standard have been defined and are being refined. New components are being added.

# Contents of this repository

This repository contains:

- A set of **example scripts** that use the PICMI standard (in the directory `Examples`)

- A **Python package** that contains the base classes for the implementation of the PICMI standard (in the directory `PICMI_Python`). In order to install this package, run:
```
cd PICMI_Python
python setup.py install
```
The latest release of this package is also available via `pip` (`pip install picmistandard``).

- The sources to generate the **Sphinx documentation** for the PICMI standard (in the directory `Docs`). In order to generate the Sphinx documentation, first install [Sphinx](), as well as the version of the PICMI standard that you would like to document (e.g. via `python setup.py`). Then type:
```
cd Docs
make html
```
You can then view the documentation by opening the file `Docs/build/html/index.html` with a standard web browser.