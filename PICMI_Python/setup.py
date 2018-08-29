#!/usr/bin/env python
# To use:
#       python setup.py install

try:
    import distutils
    from distutils.core import setup
except:
    raise SystemExit('Distutils problem')

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

setup(name = 'picmistandard',
      version = '1.0.0',
      description = 'Python base classes for PICMI standard',
      platforms = 'any',
      packages = ['picmistandard'],
      package_dir = {'picmistandard': '.'},
      cmdclass = {'build_py':build_py}
      )

