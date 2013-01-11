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

There are many functions available to calculate cosmological quantities.
They generally take a redshift as input. For example, the two cases
below give you the value of the hubble constant at z=0 (i.e., `H0`), and
the number of transverse proper kpc corresponding to an arcminute at z=3:

  >>> from astropy import cosmology
  >>> cosmology.H(0)
  70.4
  >>> cosmology.kpc_proper_per_arcmin(3)
  472.8071851564037

All the functions available are listed in the `Reference/API`_
section. These will use the "current" cosmology to calculate the
values (see `The Current Cosmology`_ section below for more
details). If you haven't set this explicitly, they will use the 9-year
WMAP cosmological parameters and print a warning message.

There are also several standard cosmologies already defined. These are
objects with methods and attributes that calculate cosmological
values. For example, the comoving distance in Mpc to redshift 4 using
the 5-year WMAP parameters:

  >>> from astropy.cosmology import WMAP5
  >>> WMAP5.comoving_distance(4)
  7329.328120760829

A full list of the pre-defined cosmologies is given by
`cosmology.parameters.available`.

An important point is that the cosmological parameters of each
instance are immutable -- that is, if you want to change, say,
`Om`, you need to make a new instance of the class.

Using `cosmology`
=================

Most of the functionality is enabled by the
`~astropy.cosmology.core.FLRW` object. This represents a
homogenous and isotropic cosmology (a cosmology characterized by the
Friedmann-Lemaitre-Robertson-Walker metric, named after the people who
solved Einstein's field equation for this special case).  However,
you can't work with this class directly, as you must specify a
dark energy model by using one of its subclasses instead,
such as `~astropy.cosmology.core.FlatLambdaCDM`.

You can create a new `~astropy.cosmology.core.FlatLambdaCDM` object with
arguments giving the hubble parameter and omega matter (both at z=0):

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  >>> cosmo
  LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

A number of additional dark energy models are provided (described below).
Note that photons and neutrinos are included in these models, so
Om0 + Ode0 is not quite one.

The pre-defined cosmologies described in the `Getting Started`_
section are instances of `~astropy.cosmology.core.FlatLambdaCDM`, and have
the same methods. So we can find the luminosity distance in Mpc to
redshift 4 by:

  >>> cosmo.luminosity_distance(4)
  35842.35374316948

or the age of the universe at z = 0 in Gyr:

  >>> cosmo.age(0)
  13.461701807287566

They also accept arrays of redshifts:

  >>> cosmo.age([0.5, 1, 1.5])
  array([ 8.42128059,  5.74698062,  4.1964541 ])

See the `~astropy.cosmology.core.FLRW` and
`~astropy.cosmology.core.FlatLambdaCDM` object docstring for all the
methods and attributes available. In addition to flat Universes,
non-flat varieties are supported such as
`~astropy.cosmology.core.LambdaCDM`.  There are also a variety of
standard cosmologies with the parameters already defined:

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.critical_density(0)       # critical density at z = 0 in g/cm^3
  9.31000313202047e-30

  >>> from astropy.cosmology import WMAP5   # WMAP 5-year
  >>> WMAP5.H(3)                    # Hubble parameter at z = 3 in km/s/Mpc
  301.71804314602889

You can see how the density parameters evolve with redshift as well

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.Om([0,1.0,2.0]), WMAP7.Ode([0.,1.0,2.0])
  (array([ 0.272     ,  0.74898525,  0.9090524 ]),
   array([ 0.72791572,  0.25055062,  0.09010261]))

Note that these don't quite add up to one even though WMAP7 assumes a
flat Universe because photons and neutrinos are included.

In addition to the `~astropy.cosmology.core.LambdaCDM` object, there
are convenience functions that calculate some of these quantities
without needing to explicitly give a cosmology - but there are more
methods available if you work directly with the cosmology object.

  >>> from astropy import cosmology
  >>> cosmology.kpc_proper_per_arcmin(3)
  472.8071851564037
  >>> cosmology.arcsec_per_kpc_proper(3)
  0.12690162477152736

These functions will perform calculations using the "current"
cosmology. This is a specific cosmology that is currently active in
`astropy` and it's described further in the following section. They
can also be explicitly given a cosmology using the `cosmo` keyword
argument. A full list of convenience functions is included below, in
the `Reference/API`_ section.


The Current Cosmology
=======================

Sometimes it's useful for Astropy functions to assume a default
cosmology so that the desired cosmology doesn't have to be specified
every time the function is called -- the convenience functions
described in the previous section are one example. For these cases
it's possible to specify a "current" cosmology.

You can set the current cosmology to a pre-defined value by using the
"default_cosmology" option in the ``[cosmology.core]`` section of the
configuration file (see :ref:`astropy_config`). Alternatively, you can
use the `~astropy.cosmology.core.set_current` function to set a
cosmology for the current Python session.

If you haven't set a current cosmology using one of the methods
described above, then the cosmology module will use the 7-year WMAP
parameters and print a warning message letting you know this. For
example, if you call a convenience function without setting the
current cosmology or using the `cosmo=` keyword you see the following
message:

  >>> from astropy import cosmology
  >>> cosmology.lookback_time(1)          # lookback time in Gyr at z=1
  WARNING: No default cosmology has been specified, using 7-year WMAP.
  [astropy.cosmology.core]
  7.787767002228743

The 9-year WMAP cosmology is also available

  >>> from astropy.cosmology import WMAP9   # WMAP 9-year
  >>> WMAP9.lookback_time(2)                # lookback time in Gyr at z=2
  10.444367272683863

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

Specifying a dark energy model
==============================

In addition to the standard `~astropy.cosmology.core.FlatLambdaCDM` model
described above, a number of additional dark energy models are
provided.  `~astropy.cosmology.core.FlatLambdaCDM` 
and `~astropy.cosmology.core.FlatLambdaCDM` assume that dark
energy is a cosmological constant, and should be the most commonly
used case.  `~astropy.cosmology.core.wCDM` assumes a constant dark
energy equation of state parameterized by :math:`w_0`. Two forms of a
variable dark energy equation of state are provided: the simple first
order linear expansion :math:`w(z) = w_0 + w_z z` by
`~astropy.cosmology.core.w0wzCDM`, as well as the common CPL form by
`~astropy.cosmology.core.w0waCDM`: :math:`w(z) = w_0 + w_a (1 - a) =
w_0 + w_a z / (1 + z)` and its generalization to include a pivot
redshift by `~astropy.cosmology.core.wpwaCDM`: :math:`w(z) = w_p + w_a
(a_p - a)`.

Users can specify their own equation of state by sub-classing
`~astropy.cosmology.core.FLRW`.  See the provided subclasses for
examples.

Relativistic Species
====================
The cosmology classes include the contribution to the energy density
from both photons and massless neutrinos.  The two parameters
controlling the proporties of these species are Tcmb0 (the temperature
of the CMB at z=0) and Neff, the effective number of neutrino species.
Both have standard default values (2.725 and 3.04, respectively; the
reason that Neff is not 3 has to do with a small bump in the neutrino
energy spectrum due to electron-positron annihilation).

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> z = [0,1.0,2.0]
  >>> WMAP7.Ogamma(z), WMAP7.Onu(z)
  (array([  4.98569503e-05,   2.74574414e-04]),
   array([  3.44204408e-05,   1.89561782e-04]),
   array([  8.42773911e-05,   4.64136197e-04]))

If you want to exclude photons and neutrinos from your calculations,
simply set the CMB Temperature to 0:

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cos = FlatLambdaCDM(70.4, 0.272, Tcmb0 = 0.0)
  >>> cos.Ogamma0, cos.Onu0
  (0.0, 0.0)

Neutrinos can be removed (while leaving photons) by setting Neff=0:

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cos = FlatLambdaCDM(70.4, 0.272, Neff=0)
  >>> cos.Ogamma([0,1,2]),cos.Onu([0,1,2])
  (array([  4.98569503e-05,   2.74623219e-04,   5.00051845e-04]),
   array([ 0.,  0.,  0.]))

While these examples used `~astropy.cosmology.core.FlatLambdaCDM`,
the above examples also apply for all of the other cosmology classes.

See Also
========

* Hogg, "Distance measures in cosmology",
  http://arxiv.org/abs/astroph/9905116
* Linder, "Exploring the Expansion History of the Universe", http://arxiv.org/abs/astro-ph/0208512
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
