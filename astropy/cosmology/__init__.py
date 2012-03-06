# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
astropy.cosmology contains routines for cosmological distance measures
and other quantities.

Capabilities will include
-------------------------

  Various cosmological densities.

  Various cosmological distance measures.

  Pre-defined sets of cosmological parameters (e.g. from WMAP).

Most of the functionality is enabled by the Cosmology object.

For example, create a new Cosmology object by given the hubble
parameter, omega matter and omega lambda (all at z=0)::

  >>> from astropy.cosmology import Cosmology
  >>> cosmo = Cosmology(H0=70, Om=0.3, Ol=0.7)
  >>> cosmo
  Cosmology(H0=70, Om=0.3, Ol=0.7, Ok=0)

Methods of the Cosmology object allow you to calculate commonly-used
quantities using your cosmology. So the comoving distance in Mpc at
redshift 4 is given by::

  >>> cosmo.comoving_distance(4)
  7170.366414463296

Age of the universe at z = 0 in Gyr::

  >>> cosmo.age(0)
  13.46698402784007

See the Cosmology object docstring for all the methods and variables
available.  When you create a new cosmology, instead of specifying H0,
Om and Ol you can use standard sets of cosmological parameters::

  # Parameters from the 7-year WMAP results
  >>> wmap7 = Cosmology('wmap7')
  # critical density at z=0 in g/cm^3
  >>> wmap7.critical_density0
  9.31000313202047e-30

  # from the 5-year results
  >>> wmap5 = Cosmology('wmap5')
  # the Hubble parameter at z = 3 in km/s/Mpc
  >>> wmap5.H(3)
  301.54148311633674

There are also several 'convenience' functions that calculate
commonly-used quantities without have to create a Cosmology
object. They use WMAP7 cosmology by default.

  >>> from astropy.cosmology import \
        kpc_proper_per_arcmin, arcsec_per_kpc_comoving
  >>> kpc_proper_per_arcmin(3)
  472.91882815884907
  >>> arcsec_per_kpc_proper(3)
  0.12687166682195736

By default these use a 7-year WMAP cosmology. The convenience
functions availableq are::

     kpc_comoving_per_arcmin
     kpc_proper_per_arcmin
     arcsec_per_kpc_comoving
     arcsec_per_kpc_proper
     distmod
     to_xyz

See their docstrings for more information.
"""
from cosmology import \
     Cosmology, kpc_comoving_per_arcmin, kpc_proper_per_arcmin, \
     arcsec_per_kpc_comoving, arcsec_per_kpc_proper, distmod, to_xyz

from parameters import WMAP7, WMAP5
