# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
astropy.cosmology contains routines for cosmological distance measures
and other quantities.

Most of the functionality in this package is enabled by the Cosmology
object. For example, create a new Cosmology object with ketwords
giving the hubble parameter, omega matter and omega lambda (all at
z=0)::

  >>> from astropy.cosmology import Cosmology
  >>> cosmo = Cosmology(H0=70, Om=0.3, Ol=0.7)
  >>> cosmo
  Cosmology(H0=70, Om=0.3, Ol=0.7, Ok=0)

The methods of this object calculate commonly used quantities with
your cosmology. So the comoving distance in Mpc at redshift 4 is given
by::

  >>> cosmo.comoving_distance(4)
  7170.366414463296

Age of the universe at z = 0 in Gyr::

  >>> cosmo.age(0)
  13.46698402784007

See the Cosmology object docstring for all the methods and variables
available.  When creating new cosmology, instead of specifying H0, Om
and Ol you can also specify a standard set of cosmological
parameters::

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

There are also several convenience functions that calculate commonly
used quantities without needing to create a Cosmology object. They use
a WMAP7 cosmology by default.

  >>> from astropy.cosmology import \
        kpc_proper_per_arcmin, arcsec_per_kpc_comoving
  >>> kpc_proper_per_arcmin(3)
  472.91882815884907
  >>> arcsec_per_kpc_proper(3)
  0.12687166682195736

By default these use a 7-year WMAP cosmology. The convenience
functions available are::

     kpc_comoving_per_arcmin
     kpc_proper_per_arcmin
     arcsec_per_kpc_comoving
     arcsec_per_kpc_proper
     distmod
     to_xyz

See their docstrings for more information.

References for most of the quantities calculated in this package are
given by Hogg (astro-ph/9905116).
"""
from cosmology import \
     Cosmology, kpc_comoving_per_arcmin, kpc_proper_per_arcmin, \
     arcsec_per_kpc_comoving, arcsec_per_kpc_proper, distmod, to_xyz
