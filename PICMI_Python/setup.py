#!/usr/bin/env python
# To use:
#       python setup.py install

from setuptools import setup

# Read the package requirements from the requirements.txt file
with open('requirements.txt') as f:
    install_requires = [ line.strip('\n') for line in f.readlines() ]

setup(name = 'picmistandard',
      version = '0.0.18',
      description = 'Python base classes for PICMI standard',
      platforms = 'any',
      packages = ['picmistandard'],
      package_dir = {'picmistandard': '.'},
      url = 'https://github.com/picmi-standard/picmi',
      install_requires = install_requires
      )
