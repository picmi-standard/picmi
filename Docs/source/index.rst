.. picmi documentation master file, created by
   sphinx-quickstart on Tue Aug 28 15:58:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Particle-In-Cell Modeling Interface (PICMI)
===========================================

Overview
--------

The Particle-In-Cell Modeling Interface (PICMI) standard establishes **conventions for the naming and structuring of input files for Particle-In-Cell simulations**.

The goal of the standard is to propose **a set (or dictionary) of names and definitions** that can be used in simulations, with as little changes as possible between codes.

Just as two persons do generally not use the exact same set of words from a dictionary (whatever the language), each code will use syntax (associated with unique functionalities of the code) that others will not. However, it is expected that, just like two persons will use and understand the same words and syntax for the most common topics, **codes will share the same language for common definitions and tasks**.

For example, it is expected that most PIC codes can share the same syntax to define a grid, a commonly-used field solver, a set of particles, etc. On the other hand, only a subset of the codes will offer, e.g. mesh refinement, field ionization, or code steering by the user. In this case, it is up to the implementation with each code to decide whether to ignore, raise a warning or an error, when an undefined statement is encountered. For more flexibility and robustness, conditional statements are available to the user to execute or read lines of the input script, based on the targeted code.

.. warning::

   The standard is still evolving at this point. The most basic components of the standard have been defined and are being refined. New components are being added.

Contents
--------

For more details on the standard and how to use it, see the links below.

.. toctree::
   :maxdepth: 1

   how_to_use/how_to_use.rst
   standard/standard.rst
