.. _astropy-cosmology:

***********************************************
Cosmological Calculations (`astropy.cosmology`)
***********************************************

.. |wCDM| replace:: :class:`~astropy.cosmology.wCDM`
.. |FlatwCDM| replace:: :class:`~astropy.cosmology.FlatwCDM`
.. |w0wzCDM| replace:: :class:`~astropy.cosmology.w0wzCDM`
.. |w0waCDM| replace:: :class:`~astropy.cosmology.w0waCDM`
.. |wpwaCDM| replace:: :class:`~astropy.cosmology.wpwaCDM`
.. |z_at_value| replace:: :func:`~astropy.cosmology.z_at_value`

Introduction
============

The :mod:`astropy.cosmology` sub-package contains classes for representing
cosmologies and utility functions for calculating commonly used quantities that
depend on a cosmological model. This includes distances, ages, and lookback
times corresponding to a measured redshift or the transverse separation
corresponding to a measured angular separation.

A number of preloaded cosmologies are available from analyses using
the WMAP and Planck satellite data. See :ref:`astropy-cosmology-realizations`.

:mod:`astropy.cosmology.units` extends the :mod:`astropy.units` sub-package,
adding and collecting cosmological units and equivalencies, like :math:`h` for
keeping track of (dimensionless) factors of the Hubble constant.

For details on reading and writing cosmologies from files, see
:ref:`cosmology_io`.

For notes on building custom Cosmology classes and interfacing
:mod:`astropy.cosmology` with 3rd-party packages, see
:ref:`astropy-cosmology-for-developers`.


Getting Started
===============

Cosmological quantities are calculated using methods of a |Cosmology| object.

Examples
--------

..
  EXAMPLE START
  Calculating Cosmological Quantities

To calculate the Hubble constant at z=0 (i.e., ``H0``) and the number of
transverse proper kiloparsecs (kpc) corresponding to an arcminute at z=3::

  >>> from astropy.cosmology import WMAP9 as cosmo
  >>> cosmo.H(0)  # doctest: +FLOAT_CMP
  <Quantity 69.32 km / (Mpc s)>

.. doctest-requires:: scipy

  >>> cosmo.kpc_proper_per_arcmin(3)  # doctest: +FLOAT_CMP
  <Quantity 472.97709620405266 kpc / arcmin>

Here |WMAP9| is a built-in object describing a cosmology with the parameters
from the nine-year WMAP results. Several other built-in cosmologies are also
available (see `Built-in Cosmologies`_). The available methods of the cosmology
object are listed in the methods summary for the |FLRW| class.

All of these methods also accept an arbitrarily-shaped array of redshifts as
input:

.. doctest-requires:: scipy

  >>> import numpy as np
  >>> from astropy.cosmology import WMAP9 as cosmo
  >>> cosmo.comoving_distance(np.array([0.5, 1.0, 1.5]))  # doctest: +FLOAT_CMP
  <Quantity [1916.06941724, 3363.07062107, 4451.7475201 ] Mpc>

You can create your own FLRW-like cosmology using one of the cosmology
classes::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
  >>> print(cosmo)  # doctest: +FLOAT_CMP
  FlatLambdaCDM(H0=70.0 km / (Mpc s), Om0=0.3, Tcmb0=2.725 K,
                Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=0.0)

Note the presence of additional cosmological parameters (e.g., ``Neff``, the
number of effective neutrino species) with default values; these can also be
specified explicitly in the call to the constructor.

..
  EXAMPLE END

The cosmology sub-package makes use of :mod:`~astropy.units`, so in many cases
returns values with units attached. Consult the documentation for that
sub-package for more details, but briefly here we will show how to access the
floating point or array values::

  >>> from astropy.cosmology import WMAP9 as cosmo
  >>> H0 = cosmo.H(0)
  >>> H0.value, H0.unit  # doctest: +FLOAT_CMP
  (np.float64(69.32), Unit("km / (Mpc s)"))


Using `astropy.cosmology`
=========================

More detailed information on using the package is provided on separate pages,
listed below.

* :ref:`astropy-cosmology-realizations`
* :ref:`Units and Equivalencies <astropy-cosmology-units-and-equivalencies>`
* :ref:`cosmology_io`
* :ref:`For Developers <astropy-cosmology-for-developers>`

Most of the functionality is enabled by the |FLRW| object. This represents a
homogeneous and isotropic cosmology (characterized by the
Friedmann-Lemaitre-Robertson-Walker metric, named after the people who solved
Einstein's field equation for this special case). However, you cannot work with
this class directly, as you must specify a dark energy model by using one of
its subclasses instead, such as |FlatLambdaCDM|.

Examples
--------

..
  EXAMPLE START
  Working with the FlatLambdaCDM Class

You can create a new |FlatLambdaCDM| object with arguments giving the Hubble
parameter and Omega matter (both at z=0)::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  >>> print(cosmo)
  FlatLambdaCDM(H0=70.0 km / (Mpc s), Om0=0.3, Tcmb0=0.0 K,
                Neff=3.04, m_nu=None, Ob0=0.0)

This can also be done more explicitly using units, which is recommended::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> import astropy.units as u
  >>> cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)


The predefined cosmologies described in the `Getting Started`_ section are
instances of |FlatLambdaCDM|, and have the same methods. So we can find the
luminosity distance to redshift 4 by:

.. doctest-requires:: scipy

  >>> cosmo.luminosity_distance(4)  # doctest: +FLOAT_CMP
  <Quantity 35842.353618623194 Mpc>

Or the age of the universe at z = 0:

.. doctest-requires:: scipy

  >>> cosmo.age(0)  # doctest: +FLOAT_CMP
  <Quantity 13.461701658024014 Gyr>

They also accept arrays of redshifts:

.. doctest-requires:: scipy

  >>> import astropy.cosmology.units as cu
  >>> cosmo.age([0.5, 1, 1.5] * cu.redshift)  # doctest: +FLOAT_CMP
  <Quantity [8.42128013, 5.74698021, 4.19645373] Gyr>

See the |FLRW| and |FlatLambdaCDM| object docstring for all of the methods and
attributes available.

..
  EXAMPLE END

..
  EXAMPLE START
  Working with Non-flat Universes with the LambdaCDM Class

In addition to flat universes, non-flat varieties are supported, such as
|LambdaCDM|. A variety of standard cosmologies with the parameters already
defined are also available (see `Built-in Cosmologies`_)
::

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.critical_density(0)  # critical density at z = 0  # doctest: +FLOAT_CMP
  <Quantity 9.31000324385361e-30 g / cm3>

You can see how the density parameters evolve with redshift as well::

  >>> import numpy as np
  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.Om(np.array([0, 1.0, 2.0]))  # doctest: +FLOAT_CMP
  array([0.272     , 0.74898522, 0.90905234])
  >>> WMAP7.Ode(np.array([0., 1.0, 2.0]))  # doctest: +FLOAT_CMP
  array([0.72791572, 0.2505506 , 0.0901026 ])

Note that these do not quite add up to one, even though |WMAP7| assumes a flat
universe, because photons and neutrinos are included. Also note that the
density parameters are unitless and so are not |Quantity| objects.

It is possible to specify the baryonic matter density at redshift zero at class
instantiation by passing the keyword argument ``Ob0``::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)
  >>> print(cosmo)
  FlatLambdaCDM(H0=70.0 km / (Mpc s), Om0=0.3, Tcmb0=0.0 K,
                Neff=3.04, m_nu=None, Ob0=0.05)

In this case the dark matter-only density at redshift 0 is available as class attribute
``Odm0`` and the redshift evolution of dark and baryonic matter densities can be
computed using the methods ``Odm`` and ``Ob``, respectively. If ``Ob0`` is not
specified, it will be set to 0 and thus ignored in further calculations.

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  >>> cosmo.Odm(1)
  np.float64(0.7741935483870968)

  >>> cosmo.Ob(1)
  np.float64(0.0)

Cosmological instances have an optional ``name`` attribute which can be used to
describe the cosmology::

  >>> from astropy.cosmology import FlatwCDM
  >>> cosmo = FlatwCDM(name='SNLS3+WMAP7', H0=71.58, Om0=0.262, w0=-1.016)
  >>> print(cosmo)
  FlatwCDM(name="SNLS3+WMAP7", H0=71.58 km / (Mpc s), Om0=0.262, Tcmb0=0.0 K, Neff=3.04,
           m_nu=None, Ob0=0.0, w0=-1.016)

..
  EXAMPLE END

This is also an example with a different model for dark energy: a flat universe
with a constant dark energy equation of state, but not necessarily a
cosmological constant. A variety of additional dark energy models are also
supported (see `Specifying a dark energy model`_).

An important point is that the cosmological parameters of each instance are
immutable — that is, if you want to change, say, ``Om``, you need to make a new
instance of the class. To make this more convenient, a
:meth:`~astropy.cosmology.Cosmology.clone` operation is provided, which allows
you to make a copy with specified values changed. Note that you cannot change
the type of cosmology with this operation (e.g., flat to non-flat).

..
  EXAMPLE START
  Making New Cosmology Instances with the .clone() Method

To make a copy of a cosmological instance using the ``clone`` operation:

  >>> from astropy.cosmology import WMAP9
  >>> newcosmo = WMAP9.clone(name='WMAP9 modified', Om0=0.3141)
  >>> WMAP9.H0, newcosmo.H0  # some values unchanged  # doctest: +FLOAT_CMP
  (<Quantity 69.32 km / (Mpc s)>, <Quantity 69.32 km / (Mpc s)>)
  >>> WMAP9.Om0, newcosmo.Om0  # some changed  # doctest: +FLOAT_CMP
  (0.2865, 0.3141)
  >>> WMAP9.Ode0, newcosmo.Ode0  # Indirectly changed since this is flat  # doctest: +FLOAT_CMP
  (np.float64(0.7134130719051658), np.float64(0.6858130719051657))

..
  EXAMPLE END

Finding the Redshift at a Given Value of a Cosmological Quantity
----------------------------------------------------------------

If you know a cosmological quantity and you want to know the redshift which it
corresponds to, you can use |z_at_value|.

Example
^^^^^^^

..
  EXAMPLE START
  Compute the Redshift at a Given Universe Age

To find the redshift using ``z_at_value``:

.. doctest-requires:: scipy

  >>> import astropy.units as u
  >>> from astropy.cosmology import Planck13, z_at_value
  >>> z_at_value(Planck13.age, 2 * u.Gyr)  # doctest: +FLOAT_CMP
  <Quantity 3.19812061 redshift>

..
  EXAMPLE END

For some quantities, there can be more than one redshift that satisfies a value.
In this case you can use the ``zmin`` and ``zmax`` keywords to restrict the
search range or set ``bracket`` to initialize it in the desired domain. See the
|z_at_value| docstring for more detailed usage examples.


Built-in Cosmologies
--------------------

A number of preloaded cosmologies are available from analyses using
the WMAP and Planck satellite data. For example:

.. doctest-requires:: scipy

  >>> from astropy.cosmology import Planck13  # Planck 2013
  >>> Planck13.lookback_time(2)  # lookback time in Gyr at z=2  # doctest: +FLOAT_CMP
  <Quantity 10.51184138 Gyr>

A full list of the predefined cosmologies can be found in
:ref:`astropy-cosmology-realizations`.


Specifying a Dark Energy Model
------------------------------

Along with the standard |FlatLambdaCDM| model described above, a number of
additional dark energy models are provided. |FlatLambdaCDM| and |LambdaCDM|
assume that dark energy is a cosmological constant, and should be the most
commonly used cases; the former assumes a flat universe, the latter allows for
spatial curvature. |FlatwCDM| and |wCDM| assume a constant dark energy equation
of state parameterized by :math:`w_{0}`. Two forms of a variable dark energy
equation of state are provided: the simple first order linear expansion
:math:`w(z) = w_{0} + w_{z} z` by |w0wzCDM|, as well as the common CPL form by
|w0waCDM|: :math:`w(z) = w_{0} + w_{a} (1 - a) = w_{0} + w_{a} z / (1 + z)`
and its generalization to include a pivot redshift by |wpwaCDM|:
:math:`w(z) = w_{p} + w_{a} (a_{p} - a)`.

Users can specify their own equation of state by subclassing |FLRW|. See the
provided subclasses for examples. It is advisable to stick to subclassing
|FLRW| rather than one of its subclasses, since some of them use internal
optimizations that also need to be propagated to any  subclasses. Users wishing
to use similar tricks (which can make distance calculations much faster) should
consult the cosmology module source code for details.

Photons and Neutrinos
---------------------

The cosmology classes (can) include the contribution to the energy density from
both photons and neutrinos. By default, the latter are assumed massless. The
three parameters controlling the properties of these species, which are
arguments to the initializers of all of the cosmological classes, are ``Tcmb0``
(the temperature of the cosmic microwave background at z=0), ``Neff`` (the
effective number of neutrino species), and ``m_nu`` (the rest mass of the
neutrino species). ``Tcmb0`` and ``m_nu`` should be expressed as unit
Quantities. All three have standard default values — 0 K, 3.04, and 0 eV,
respectively. (The reason that ``Neff`` is not 3 has to do primarily with a
small bump in the neutrino energy spectrum due to electron-positron
annihilation, but is also affected by weak interaction physics.) Setting the
CMB temperature to 0 removes the contribution of both neutrinos and photons.
This is the default to ensure these components are excluded unless the user
explicitly requests them.

Massive neutrinos are treated using the approach described in the
WMAP seven-year cosmology paper (Komatsu et al. 2011, ApJS, 192, 18, section
3.3). This is not the simple
:math:`\Omega_{\nu 0} h^2 = \sum_i m_{\nu\, i} / 93.04\,\mathrm{eV}`
approximation. Also note that the values of :math:`\Omega_{\nu}(z)` include
both the kinetic energy and the rest mass energy components, and that the
|Planck13| and |Planck15| cosmologies include a single species of neutrinos
with non-zero mass (which is not included in :math:`\Omega_{m0}`).

Adding massive neutrinos can have significant performance implications. In
particular, the computation of distance measures and lookback times are factors
of three to four times slower than in the massless neutrino case. Therefore, if
you need to compute many distances in such a cosmology and performance is
critical, it is particularly useful to calculate them on a grid and use
interpolation.

Examples
^^^^^^^^

..
  EXAMPLE START
  Calculating the Contribution of Photons and Neutrinos to the Energy Density

The contribution of photons and neutrinos to the total mass-energy density can
be found as a function of redshift::

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.Ogamma0, WMAP7.Onu0  # Current epoch values  # doctest: +FLOAT_CMP
  (np.float64(4.985694972799396e-05), np.float64(3.442154948307989e-05))
  >>> z = np.array([0, 1.0, 2.0])
  >>> WMAP7.Ogamma(z), WMAP7.Onu(z)  # doctest: +FLOAT_CMP
  (array([4.98603986e-05, 2.74593395e-04, 4.99915942e-04]),
   array([3.44239306e-05, 1.89580995e-04, 3.45145089e-04]))

If you want to exclude photons and neutrinos from your calculations, you can
set ``Tcmb0`` to 0 (which is also the default)::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> import astropy.units as u
  >>> cos = FlatLambdaCDM(70.4 * u.km / u.s / u.Mpc, 0.272, Tcmb0 = 0.0 * u.K)
  >>> cos.Ogamma0, cos.Onu0
  (np.float64(0.0), np.float64(0.0))

You can include photons but exclude any contributions from neutrinos by setting
``Tcmb0`` to be non-zero (2.725 K is the standard value for our Universe) but
setting ``Neff`` to 0::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cos = FlatLambdaCDM(70.4, 0.272, Tcmb0=2.725, Neff=0)
  >>> cos.Ogamma(np.array([0, 1, 2]))  # Photons are still present  # doctest: +FLOAT_CMP
  array([4.98603986e-05, 2.74642208e-04, 5.00086413e-04])
  >>> cos.Onu(np.array([0, 1, 2]))  # But not neutrinos  # doctest: +FLOAT_CMP
  array([0., 0., 0.])

The number of neutrino species is assumed to be the floor of ``Neff``, which in
the default case is ``Neff=3``. Therefore, if non-zero neutrino masses are
desired, then three masses should be provided. However, if only one value is
provided, all of the species are assumed to have the same mass. ``Neff`` is
assumed to be shared equally between each species.

::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> import astropy.units as u
  >>> H0 = 70.4 * u.km / u.s / u.Mpc
  >>> m_nu = 0 * u.eV
  >>> cosmo = FlatLambdaCDM(H0, 0.272, Tcmb0=2.725, m_nu=m_nu)
  >>> cosmo.has_massive_nu
  False
  >>> cosmo.m_nu  # doctest: +FLOAT_CMP
  <Quantity [0., 0., 0.] eV>
  >>> m_nu = [0.0, 0.05, 0.10] * u.eV
  >>> cosmo = FlatLambdaCDM(H0, 0.272, Tcmb0=2.725, m_nu=m_nu)
  >>> cosmo.has_massive_nu
  True
  >>> cosmo.m_nu  # doctest: +FLOAT_CMP
  <Quantity [0.  , 0.05, 0.1 ] eV>
  >>> cosmo.Onu(np.array([0, 1.0, 15.0]))  # doctest: +FLOAT_CMP
  array([0.00327011, 0.00896845, 0.01257946])
  >>> cosmo.Onu(1) * cosmo.critical_density(1)  # doctest: +FLOAT_CMP
  <Quantity 2.444380380370406e-31 g / cm3>

While these examples used |FlatLambdaCDM|, the above examples also apply for
all of the other cosmology classes.

..
  EXAMPLE END

See Also
========

* Hogg, "Distance measures in cosmology",
  https://arxiv.org/abs/astro-ph/9905116
* Linder, "Exploring the Expansion History of the Universe", https://arxiv.org/abs/astro-ph/0208512
* NASA's Legacy Archive for Microwave Background Data Analysis,
  https://lambda.gsfc.nasa.gov/

Range of Validity and Reliability
=================================

The code in this sub-package is tested against several widely used online
cosmology calculators and has been used to perform many calculations in
refereed papers. You can check the range of redshifts over which the code is
regularly tested in the test suite. If you find any bugs, please let us know
by `opening an issue at the GitHub repository
<https://github.com/astropy/astropy/issues>`_!

A more difficult question is the range of redshifts over which the code is
expected to return valid results. This is necessarily model-dependent, but in
general you should not expect the numeric results to be well behaved for
redshifts more than a few times larger than the epoch of matter-radiation
equality (so, for typical models, not above z = 5-6,000, but for some models
much lower redshifts may be ill-behaved). In particular, you should pay
attention to warnings from the :mod:`scipy.integrate` package about integrals
failing to converge (which may only be issued once per session).

The built-in cosmologies use the parameters as listed in the respective papers.
These provide only a limited range of precision, and so you should not expect
derived quantities to match beyond that precision. For example, the Planck 2013
and 2015 results only provide the Hubble constant to four digits. Therefore,
they should not be expected to match the age quoted by the Planck team to
better than that, despite the fact that five digits are quoted in the papers.

Reference/API
=============

.. toctree::
   :maxdepth: 2

   realizations
   Units and Equivalencies <units>
   io/index
   For Developers <dev>
   ref_api
