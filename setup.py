#!/usr/bin/python
from distutils.core import setup, Extension
pySpinModule = Extension('pySpin', sources = ['pySpin.c','random.c'])
setup (name = 'PyKMC',
       version = '0.1',
       description = "This is a package for doing Kinetic Monte Carlo on Lattice Systems",
       ext_modules = [pySpinModule] )
