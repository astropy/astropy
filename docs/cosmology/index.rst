***********************************************
Cosmological Calculations (`astropy.cosmology`)
***********************************************

Introduction
============

The `astropy.cosmology` subpackage contains classes for representing
cosmologies, and utility functions for calculating commonly used
quantities that depend on a cosmological model. This includes
distances, ages and lookback times corresponding to a measured
redshift or the transverse separation corresponding to a measured
angular separation.


Getting Started
===============

Create a `FLRWCosmology` object with arguments giving the hubble
parameter, omega matter and omega lambda (all at z=0):

  >>> from astropy.cosmology import FLRWCosmology
  >>> cosmo = FLRWCosmology(H0=70, Om=0.3, Ol=0.7)
  >>> cosmo
  FLRWCosmology(H0=70, Om=0.3, Ol=0.7, Ok=0)

Methods of this object calculate commonly used quantities with your
cosmology. For example, the comoving distance in Mpc at redshift 4 is
given by:

  >>> cosmo.comoving_distance(4)
  7170.366414463296

These methods also accept arrays of redshifts:

  >>> cosmo.comoving_distance([4.1, 4.2, 4.3])
  array([ 7238.65549685,  7304.99885736,  7369.48598852])

There are several standard cosmologies already defined:

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.comoving_distance(4)
  7352.203452009956

A full list of the pre-defined cosmologies is given by
`cosmology.parameters.available`.

There are also convenience functions that calculate quantities without
needing to create a Cosmology object.

  >>> from astropy import cosmology
  >>> cosmology.kpc_proper_per_arcmin(3)
  472.91882815884907



Using `cosmology`
=================

Most of the functionality is enabled by the
`~astropy.cosmology.core.FLRWCosmology` object. This represents a
homogenous and isotropic cosmology (a cosmology characterized by the
Friedmann-Lemaitre-Robertson-Walker metric, named after the people who
solved Einstein's field equation for this special case).

You can create a new `FLRWCosmology` object with arguments giving the
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
variables available. There are also a variety of standard cosmologies
with the parameters already defined:

  >>> from cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.critical_density(0)       # critical density at z = 0 in g/cm^3
  9.31000313202047e-30

  >>> from cosmology import WMAP5   # WMAP 5-year
  >>> WMAP5.H(3)                    # Hubble parameter at z = 3 in km/s/Mpc
  301.54148311633674


In addition to the `FLRWCosmology` object, there are convenience
functions that calculate quantities without needing to explicitly give
a cosmology.

  >>> from astropy import cosmology
  >>> cosmology.kpc_proper_per_arcmin(3)
  472.91882815884907
  >>> cosmology.arcsec_per_kpc_proper(3)
  0.12687166682195736

These functions will perform calculations using the "current"
cosmology. This is a specific cosmology that is currently active in
`astropy` and it's described further in the following section. They
can also be explicitly given a cosmology using the `cosmo` keyword
argument. A full list of convenience functions is included below, in
the `Reference/API`_ section.


The "Current" Cosmology
=======================

Sometimes it's useful for Astropy functions to assume a default
cosmology so that the desired cosmology doesn't have to be specified
every time the function is called.  The convenience functions from the
previous section are one such example. For these cases it's possible
to specify a "current" cosmology.

You can set the current cosmology to a pre-defined value by using the
"default_cosmology" option in the ``[cosmology.core]`` section of the
configuration file (see :ref:`astropy_config`). Alternatively, you can
use the `~astropy.cosmology.core.set_current` function to set a
cosmology for the current Python session.

If you haven't set a current cosmology using one of the methods
described above, then the cosmology module will use `WMAP7` and print
a warning message letting you know this. For example, if you call a
convenience function without setting the current cosmology or using
the `cosmo=` keyword you see the following message:

  >>> from astropy import cosmology
  >>> cosmology.lookback_time(1)          # lookback time in Gyr at z=1
  WARNING: No default cosmology has been specified, using 7-year WMAP. 
  [astropy.cosmology.core]
  7.788414051773566

.. note::

    In general it's better to use an explicit cosmology (for example
    ``WMAP7.H(0)`` instead of ``cosmology.H(0)``). The motivation for
    this is that when you go back to use the code at a later date or
    share your scripts with someone else, the default cosmology may
    have changed. Use of the convenience functions should generally be
    reserved for interactive work or cases where the flexibility of
    quickly changing between different cosmologies is for some reason
    useful. Alternatively, putting (for example)
    ``cosmology.set_current(WMAP7)`` at the top of your code will
    ensure that the right cosmology is always used.


Using `cosmology` inside Astropy
--------------------------------

If you are writing code for the `astropy` core or an affiliated package,
it is strongly recommended that you use the the current cosmology
through the `~astropy.cosmology.core.get_current` function. It is also
recommended that you provide an override option something like the
following::

    def myfunc(..., cosmo=None):
	from astropy.cosmology import get_current

	if cosmo is None:
	    cosmo = get_current()

	... your code here ...

This ensures that all code consistently uses the current cosmology
unless explicitly overridden.


See Also
========

* Hogg, "Distance measures in cosmology",
  http://arxiv.org/abs/astroph/9905116
* NASA's Legacy Archive for Microwave Background Data Analysis,
  http://lambda.gsfc.nasa.gov/

Range of validity and reliability
=================================

The code in this sub-package is tested against several widely-used
online cosmology calculators, and has been used to perform
calculations in refereed papers. You can check the range of redshifts
over which the code is regularly tested in the module
`astropy.cosmology.tests.test_cosmology`. Note that the energy density
due to radiation is assumed to be negligible, which is valid for
redshifts less than about 10. If you find any bugs, please let us know
by `opening an issue at the github repository
<https://github.com/astropy/astropy/issues>`_!

Reference/API
=============

.. automodapi:: astropy.cosmology
