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


# For developers

This directory is the Python version of the PICMI standard. The files contain the classes that form the standard. These files
consist of picmistandard package.

To submit changes, create a fork of the repo on github and create a new branch in the fork where the changes are to be made.
When the changes are ready submit a pull request on github for review. If acceptable, the request will be merged in.

Some packages specify the version of the standard that is to be used. For new features to be included, the version of picmi needs to
be updated. The update has several steps, updating the version number, adding a tag, pushing the changes, and updating the version
of PyPI.

Currently, the version number is modified by hand by modifying the number in the setup.py file, incrementing the third
number. After changing the file, make a commit with the change with a comment like "Version 0.0.16", updating the number of course.
Please don't make any other changes in the commit. Note that it is ok to push a version update directly to the main branch assuming
it has been approved.

Next, add a tag and push the changes to the repo:

  ```
  git tag -a 0.0.16 -m "Release 0.0.16"
  git push origin
  git push origin --tags
  ```

Change the tag number as needed.

The final step is to update the version on PyPI. It is recommended to use the twine command. Here are the commands:

  ```
  python setup.py sdist bdist_wheel
  twine upload dist/*
  ```

The first builds the distribution, the second uploads it. Note that it is good practice to delete the dist directory
afterwards so that each time uploads only the most recent version.

