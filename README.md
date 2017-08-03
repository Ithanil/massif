# massif
Modules for Atomistic Simulation Software in Fortran

This project is meant to be a simple collection of Modern Fortran modules that provide basic functionality for various kinds of atomistic many body system simulations (e.g. Molecular Dynamics or Metropolis Monte Carlo).

Currently it contains the following modules:

kinds: Stores the size of single and double precision floats (+complex type).

natconst: Stores various frequently used natural constants.

mcwalkers: Object-Oriented implementation of walker for Monte Carlo (MC) algorithms.

mcinterfaces: Provides an interface for sampling of observables via MC.

An old file concerning PIMD is also included, but needs to be revisited:

sys: Provides a data structure to represent the coordinate state of molecular system in the case of PIMD.
