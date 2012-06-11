Cosmological Calculations (`astropy.packagename`)
=================================================

Introduction
------------

The `astropy.cosmology` package contains classes for representing cosmologies,
and utility functions for common cosmological calculations.


Getting Started
---------------

This section needs to be populated.


Using `cosmology`
-----------------

The cosmology package allows you to calculate many commonly used
quantities that depend on the cosmological model, such as distances,
ages and lookback times corresponding to a measured redshift or the
transverse separation corresponding to a measure angular separation.

Most of the functionality is enabled by the `FLRWCosmology`
object. This represents a homgeneous and isotropic cosmology (also
known as a Friedmann-Lemaitre-Robertson-Walker cosmology after the
people who solved Einstein's field equation for this special case).

To create a new `FLRWCosmology` object with arguments giving the
hubble parameter, omega matter and omega lambda (all at z=0):

  >>> from astropy.cosmology import FLRWCosmology
  >>> cosmo = FLRWCosmology(H0=70, Om=0.3, Ol=0.7)
  >>> cosmo
  FLRWCosmology(H0=70, Om=0.3, Ol=0.7, Ok=0)

The methods of this object calculate commonly used quantities with
your cosmology. For example, the comoving distance in Mpc at redshift
4 is given by:

  >>> cosmo.comoving_distance(4)
  7170.366414463296

The age of the universe at z = 0 in Gyr:

  >>> cosmo.age(0)
  13.46698402784007

See the `FLRWCosmology` object docstring for all the methods and
variables available.  There are several standard cosmologies already
defined:

  >>> from cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.critical_density0       # critical density at z = 0 in g/cm^3
  9.31000313202047e-30

  >>> from cosmology import WMAP5   # WMAP 5-year
  >>> WMAP5.H(3)                    # Hubble parameter at z = 3 in km/s/Mpc
  301.54148311633674

There is also a 'current' cosmology that will be used by relevant
astropy functions if no cosmology instance is explicitly passed to
them. This can be set with `set_current()` or a configuration file
option (option name 'default_cosmology'), and is accessed with
`get_current()`. If you don't set the current explicitly, either with
`set_current()` or with the configuration file option, then
`get_current()` returns the 7 year WMAP7 cosmology and a warning
message is printed.

There are also several convenience functions that calculate quantities
without needing to create a Cosmology object.

  >>> from astropy import cosmology
  >>> cosmology.kpc_proper_per_arcmin(3)
  472.91882815884907
  >>> cosmology.arcsec_per_kpc_proper(3)
  0.12687166682195736

These use the current cosmology, unless overridden by a `cosmo=`
keyword argument. The convenience functions available are:

     `~astropy.cosmology.kpc_comoving_per_arcmin`
     `~astropy.cosmology.kpc_proper_per_arcmin`
     `~astropy.cosmology.arcsec_per_kpc_comoving`
     `~astropy.cosmology.arcsec_per_kpc_proper`
     `~astropy.cosmology.distmod`

See their docstrings for more information.

References for most of the quantities calculated in this package are
given by Hogg (astro-ph/9905116).


See Also
--------

* Hogg, "Distance measures in cosmology", http://arxiv.org/abs/astroph/9905116
* NASA's Legacy Archive for Microwave Background Data Analysis, http://lambda.gsfc.nasa.gov/

Reference/API
-------------

.. automodapi:: astropy.cosmology


Acknowledgments and Licenses (optional)
---------------------------------------

Any acknowledgements or licenses needed for this package - remove the
section if none are necessary.
