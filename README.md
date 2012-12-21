PyKMC
=====

Kinetic Monte Carlo code, specifically targeted to simulating lattice glass models

Info
-------
Started Monday, October 8, 2012 by Glen Hocky

Authors
-------
Glen Hocky

Quick Instructions
------------------
To build:
    python setup.py build_ext -i

To run:
 * numpy should be in PYTHONPATH
 * add this directory to PYTHONPATH

Contents
-----------
* README.md

    This file, describing the overall contents of this code package
* fileio/

    This should contain utilities for writing and reading configurations and
        trajectories
* models/

    This should contain rules for the dynamics, and code that builds the list
        of potential moves for the current configuration
    This now contains lattice definitions and tools for determining neighbors
* tests/

    This should contain scripts that run the code in known configurations to 
        check for expected results

