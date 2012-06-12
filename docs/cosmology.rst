Cosmological Calculations (`astropy.cosmology`)
===============================================

Introduction
------------

The `astropy.cosmology` subpackage contains classes for representing
cosmologies, and utility functions for calculate many commonly used
quantities that depend on the cosmological model. This includes
distances, ages and lookback times corresponding to a measured redshift
or the transverse separation corresponding to a measure angular
separation.

An important concept in `astropy.cosmology` is the "current" cosmology.
This is the specific cosmological model and choice of parameters that are
currently active in `astropy`. Other parts of Astropy that require
knowledge of the background cosmology will use this cosmology in their
calculations to maintain consistency. See `Getting Started`_ for a
description of how to change the current cosmology that is in use.




Getting Started
---------------

To do a calculation defined in one of the convinience functions, you can
simply call the function with the relevant redshift::

    >>> from astropy import cosmology
    >>> cosmology.distmod(0.5)
    WARNING: No default cosmology has been specified, using 7-year WMAP. [astropy.cosmology.core]
    42.270203330485998

Note that calling these functions without specifying a cosmology will
cause a warning to appear and the default (WMAP7) will be adopted. You
can get rid of this by specifying a default cosmology or setting the
current cosmology directly. The default current cosmology can be changed
by changing the "default_cosmology" option in the ``[cosmology.core]``
section of the configuration file to your preferred cosmology (see
:ref:`astropy_config`). Alternatively, you can use the
`~astropy.cosmology.core.set_current`. function to specify a cosmology
for use in the current python session.

Most of the other functionality is implmented as either methods or
attributes of the current cosmology object. Use
`~astropy.cosmology.core.get_current` to get this object::

    >>> from astropy.cosmology import get_current
    >>> cosmo = get_current()
    >>> cosmo.h
    0.704
    >>> cosmo.lookback_time(1)
    7.788414051773566
    >>> cosmo.critical_density(0)
    9.3100031320204701e-30
    >>> cosmo.critical_density(0.5)
    1.5324265155305696e-29



Using `cosmology`
-----------------

Most of the functionality is enabled by the
`~astropy.cosmology.core.FLRWCosmology` object. This represents a
homgeneous and isotropic cosmology (a cosmology characterized by the
Friedmann-Lemaitre-Robertson-Walker metric after the people who solved
Einstein's field equation for this special case).

While `astropy.cosmology` includes a variety of standard cosmologies
with the parameters already defined (see below), you can create a new
`FLRWCosmology` object with arguments giving the hubble parameter, omega
matter and omega lambda (all at z=0):

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


There are also several convenience functions that calculate quantities
without needing to create a Cosmology object.

  >>> from astropy import cosmology
  >>> cosmology.kpc_proper_per_arcmin(3)
  472.91882815884907
  >>> cosmology.arcsec_per_kpc_proper(3)
  0.12687166682195736

These use the current cosmology, unless overridden by a `cosmo=` keyword
argument. A full list of convinience functions is included below, in the
`Reference/API`_ section.


Using `cosmology` inside Astropy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are writing code for the `astropy` core or an affiliated package,
it is strongly recommended that you use the the current cosmology
through the `~astropy.cosmology.core.get_current` function. It also also
recommended that you provide an override option something like the
following::

    def myfunc(...,cosmo=None):
        from astropy.cosmology import get_current

        if cosmo is None:
            cosmo = get_current()

        ... your code here ...

This ensures that all code consistently uses the current cosmology
unless explicitly overridden.


See Also
--------

* Hogg, "Distance measures in cosmology", http://arxiv.org/abs/astroph/9905116
* NASA's Legacy Archive for Microwave Background Data Analysis, http://lambda.gsfc.nasa.gov/

Reference/API
-------------

.. automodapi:: astropy.cosmology

