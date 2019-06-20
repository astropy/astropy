.. _astropy-cosmology:

***********************************************
Cosmological Calculations (`astropy.cosmology`)
***********************************************

Introduction
============

The `astropy.cosmology` sub-package contains classes for representing
cosmologies and utility functions for calculating commonly used
quantities that depend on a cosmological model. This includes
distances, ages, and lookback times corresponding to a measured
redshift or the transverse separation corresponding to a measured
angular separation.


Getting Started
===============

Cosmological quantities are calculated using methods of a
:class:`~astropy.cosmology.Cosmology` object. For example, to calculate the
Hubble constant at z=0 (i.e., ``H0``) and the number of transverse proper
kiloparsecs (kpc) corresponding to an arcminute at z=3::

  >>> from astropy.cosmology import WMAP9 as cosmo
  >>> cosmo.H(0)  # doctest: +FLOAT_CMP
  <Quantity 69.32 km / (Mpc s)>

.. doctest-requires:: scipy

  >>> cosmo.kpc_proper_per_arcmin(3)  # doctest: +FLOAT_CMP
  <Quantity 472.97709620405266 kpc / arcmin>

Here WMAP9 is a built-in object describing a cosmology with the
parameters from the nine-year WMAP results. Several other built-in
cosmologies are also available (see `Built-in Cosmologies`_). The
available methods of the cosmology object are listed in the methods
summary for the `~astropy.cosmology.FLRW` class. If you are using
IPython you can also use tab completion to print a list of the
available methods. To do this, after importing the cosmology as in the
above example, type ``cosmo.`` at the IPython prompt and then press
the tab key.

All of these methods also accept an arbitrarily-shaped array of
redshifts as input:

.. doctest-requires:: scipy

  >>> from astropy.cosmology import WMAP9 as cosmo
  >>> cosmo.comoving_distance([0.5, 1.0, 1.5])  # doctest: +FLOAT_CMP
  <Quantity [1916.06941724, 3363.07062107, 4451.7475201 ] Mpc>

You can create your own FLRW-like cosmology using one of the cosmology
classes::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
  >>> cosmo  # doctest: +FLOAT_CMP
  FlatLambdaCDM(H0=70 km / (Mpc s), Om0=0.3, Tcmb0=2.725 K,
                Neff=3.04, m_nu=[0. 0. 0.] eV, Ob0=None)

Note the presence of additional cosmological parameters (e.g., ``Neff``,
the number of effective neutrino species) with default values; these
can also be specified explicitly in the call to the constructor.

The cosmology sub-package makes use of `~astropy.units`, so in many
cases returns values with units attached. Consult the documentation
for that sub-package for more details, but briefly here we will show how to
access the floating point or array values::

  >>> from astropy.cosmology import WMAP9 as cosmo
  >>> H0 = cosmo.H(0)
  >>> H0.value, H0.unit  # doctest: +FLOAT_CMP
  (69.32, Unit("km / (Mpc s)"))


Using `astropy.cosmology`
=========================

Most of the functionality is enabled by the `~astropy.cosmology.FLRW`
object. This represents a homogeneous and isotropic cosmology
(characterized by the Friedmann-Lemaitre-Robertson-Walker metric,
named after the people who solved Einstein's field equation for this
special case). However, you cannot work with this class directly, as
you must specify a dark energy model by using one of its subclasses
instead, such as `~astropy.cosmology.FlatLambdaCDM`.

You can create a new `~astropy.cosmology.FlatLambdaCDM` object with
arguments giving the Hubble parameter and Omega matter (both at z=0)::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  >>> cosmo
  FlatLambdaCDM(H0=70 km / (Mpc s), Om0=0.3, Tcmb0=0 K,
                Neff=3.04, m_nu=None, Ob0=None)

This can also be done more explicitly using units, which is recommended::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> import astropy.units as u
  >>> cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

However, most of the parameters that accept units (``H0``, ``Tcmb0``)
have default units, so unit quantities do not have to be used (with the
exception of neutrino masses, where you must supply a unit if you want massive
neutrinos).

The predefined cosmologies described in the `Getting Started`_
section are instances of `~astropy.cosmology.FlatLambdaCDM`, and have
the same methods. So we can find the luminosity distance to
redshift 4 by:

.. doctest-requires:: scipy

  >>> cosmo.luminosity_distance(4)  # doctest: +FLOAT_CMP
  <Quantity 35842.353618623194 Mpc>

Or the age of the universe at z = 0:

.. doctest-requires:: scipy

  >>> cosmo.age(0)  # doctest: +FLOAT_CMP
  <Quantity 13.461701658024014 Gyr>

They also accept arrays of redshifts:

.. doctest-requires:: scipy

  >>> cosmo.age([0.5, 1, 1.5]).value  # doctest: +FLOAT_CMP
  array([8.42128013, 5.74698021, 4.19645373])

See the `~astropy.cosmology.FLRW` and
`~astropy.cosmology.FlatLambdaCDM` object docstring for all of the
methods and attributes available. In addition to flat universes,
non-flat varieties are supported, such as
`~astropy.cosmology.LambdaCDM`. A variety of standard
cosmologies with the parameters already defined are also available
(see `Built-in Cosmologies`_)::

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.critical_density(0)  # critical density at z = 0  # doctest: +FLOAT_CMP
  <Quantity 9.31000324385361e-30 g / cm3>

You can see how the density parameters evolve with redshift as well::

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.Om([0, 1.0, 2.0]), WMAP7.Ode([0., 1.0, 2.0])  # doctest: +FLOAT_CMP
  (array([0.272     , 0.74898522, 0.90905234]),
   array([0.72791572, 0.2505506 , 0.0901026 ]))

Note that these do not quite add up to one, even though WMAP7 assumes a
flat universe, because photons and neutrinos are included. Also note
that the density parameters are unitless and so are not
`~astropy.units.Quantity` objects.

It is possible to specify the baryonic matter density at redshift zero
at class instantiation by passing the keyword argument ``Ob0``::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)
  >>> cosmo
  FlatLambdaCDM(H0=70 km / (Mpc s), Om0=0.3, Tcmb0=0 K,
                Neff=3.04, m_nu=None, Ob0=0.05)

In this case the dark matter-only density at redshift 0 is
available as class attribute ``Odm0`` and the redshift evolution of
dark and baryonic matter densities can be computed using the methods
``Odm`` and ``Ob``, respectively. If ``Ob0`` is not specified at class
instantiation, it defaults to ``None`` and any method relying on it
being specified will raise a ``ValueError``:

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  >>> cosmo.Odm(1)
  Traceback (most recent call last):
  ...
  ValueError: Baryonic density not set for this cosmology, unclear
  meaning of dark matter density

Cosmological instances have an optional ``name`` attribute which can be
used to describe the cosmology::

  >>> from astropy.cosmology import FlatwCDM
  >>> cosmo = FlatwCDM(name='SNLS3+WMAP7', H0=71.58, Om0=0.262, w0=-1.016)
  >>> cosmo
  FlatwCDM(name="SNLS3+WMAP7", H0=71.6 km / (Mpc s), Om0=0.262,
           w0=-1.02, Tcmb0=0 K, Neff=3.04, m_nu=None, Ob0=None)

This is also an example with a different model for dark energy: a flat
universe with a constant dark energy equation of state, but not
necessarily a cosmological constant. A variety of additional dark
energy models are also supported (see `Specifying a dark energy
model`_).

An important point is that the cosmological parameters of each
instance are immutable — that is, if you want to change, say,
``Om``, you need to make a new instance of the class. To make
this more convenient, a ``clone`` operation is provided, which
allows you to make a copy with specified values changed.
Note that you cannot change the type of cosmology with this operation
(e.g., flat to non-flat). For example:

  >>> from astropy.cosmology import WMAP9
  >>> newcosmo = WMAP9.clone(name='WMAP9 modified', Om0=0.3141)
  >>> WMAP9.H0, newcosmo.H0  # some values unchanged  # doctest: +FLOAT_CMP
  (<Quantity 69.32 km / (Mpc s)>, <Quantity 69.32 km / (Mpc s)>)
  >>> WMAP9.Om0, newcosmo.Om0  # some changed  # doctest: +FLOAT_CMP
  (0.2865, 0.3141)
  >>> WMAP9.Ode0, newcosmo.Ode0  # Indirectly changed since this is flat  # doctest: +FLOAT_CMP
  (0.7134130719051658, 0.6858130719051657)

Finding the Redshift at a Given Value of a Cosmological Quantity
----------------------------------------------------------------

If you know a cosmological quantity and you want to know the
redshift which it corresponds to, you can use ``z_at_value``:

.. doctest-requires:: scipy

  >>> import astropy.units as u
  >>> from astropy.cosmology import Planck13, z_at_value
  >>> z_at_value(Planck13.age, 2 * u.Gyr)  # doctest: +FLOAT_CMP
  3.1981226843560968

For some quantities, there can be more than one redshift that satisfies
a value. In this case you can use the ``zmin`` and ``zmax`` keywords
to restrict the search range. See the ``z_at_value`` docstring for more
detailed usage examples.


Built-in Cosmologies
--------------------

A number of preloaded cosmologies are available from analyses using
the WMAP and Planck satellite data. For example:

.. doctest-requires:: scipy

  >>> from astropy.cosmology import Planck13  # Planck 2013
  >>> Planck13.lookback_time(2)  # lookback time in Gyr at z=2  # doctest: +FLOAT_CMP
  <Quantity 10.51184138 Gyr>

A full list of the predefined cosmologies is given by
``cosmology.parameters.available`` and summarized below:

========  ============================== ====  ===== =======
Name      Source                         H0    Om    Flat
========  ============================== ====  ===== =======
WMAP5     Komatsu et al. 2009            70.2  0.277 Yes
WMAP7     Komatsu et al. 2011            70.4  0.272 Yes
WMAP9     Hinshaw et al. 2013            69.3  0.287 Yes
Planck13  Planck Collab 2013, Paper XVI  67.8  0.307 Yes
Planck15  Planck Collab 2015, Paper XIII 67.7  0.307 Yes
========  ============================== ====  ===== =======

Currently, all are instances of `~astropy.cosmology.FlatLambdaCDM`.
More details about exactly where each set of parameters comes from
are available in the docstring for each object::

  >>> from astropy.cosmology import WMAP7
  >>> print(WMAP7.__doc__)
  WMAP7 instance of FlatLambdaCDM cosmology
  (from Komatsu et al. 2011, ApJS, 192, 18, doi: 10.1088/0067-0049/192/2/18.
  Table 1 (WMAP + BAO + H0 ML).)


Specifying a Dark Energy Model
------------------------------

Along with the standard `~astropy.cosmology.FlatLambdaCDM` model
described above, a number of additional dark energy models are
provided. `~astropy.cosmology.FlatLambdaCDM`
and `~astropy.cosmology.LambdaCDM` assume that dark
energy is a cosmological constant, and should be the most commonly
used cases; the former assumes a flat universe, the latter allows
for spatial curvature. `~astropy.cosmology.FlatwCDM` and
`~astropy.cosmology.wCDM` assume a constant dark
energy equation of state parameterized by :math:`w_{0}`.
Two forms of a variable dark energy equation of state are provided: the simple
first order linear expansion :math:`w(z) = w_{0} + w_{z} z` by
`~astropy.cosmology.w0wzCDM`, as well as the common CPL form by
`~astropy.cosmology.w0waCDM`: :math:`w(z) = w_{0} + w_{a} (1 - a) =
w_{0} + w_{a} z / (1 + z)` and its generalization to include a pivot
redshift by `~astropy.cosmology.wpwaCDM`: :math:`w(z) = w_{p} + w_{a}
(a_{p} - a)`.

Users can specify their own equation of state by subclassing
`~astropy.cosmology.FLRW`. See the provided subclasses for
examples. It is recommended, but not required, that all arguments to the
constructor of a new subclass be available as properties, since the
``clone`` method assumes this is the case. It is also advisable
to stick to subclassing `~astropy.cosmology.FLRW` rather than one of
its subclasses, since some of them use internal optimizations that
also need to be propagated to any subclasses. Users wishing to
use similar tricks (which can make distance calculations much faster)
should consult the cosmology module source code for details.

Photons and Neutrinos
---------------------

The cosmology classes (can) include the contribution to the energy density
from both photons and neutrinos. By default, the latter are assumed
massless. The three parameters controlling the properties of these
species, which are arguments to the initializers of all of the
cosmological classes, are ``Tcmb0`` (the temperature of the cosmic microwave
background at z=0), ``Neff`` (the effective number of neutrino species), and
``m_nu`` (the rest mass of the neutrino species). ``Tcmb0`` and ``m_nu`` should
be expressed as unit Quantities.
All three have standard default values — 0 K, 3.04, and 0 eV, respectively.
(The reason that ``Neff`` is not 3 has to do primarily with a small bump in the
neutrino energy spectrum due to electron-positron annihilation, but is also
affected by weak interaction physics.) Setting the CMB temperature to 0
removes the contribution of both neutrinos and photons. This is the default to
ensure these components are excluded unless the user explicitly requests them.

Massive neutrinos are treated using the approach described in the
WMAP seven-year cosmology paper (Komatsu et al. 2011, ApJS, 192, 18, section
3.3). This is not the simple
:math:`\Omega_{\nu 0} h^2 = \sum_i m_{\nu\, i} / 93.04\,\mathrm{eV}`
approximation. Also note that the values of :math:`\Omega_{\nu}(z)`
include both the kinetic energy and the rest mass energy components,
and that the Planck13 and Planck15 cosmologies include a single
species of neutrinos with non-zero mass (which is not included in
:math:`\Omega_{m0}`).

Adding massive neutrinos can have significant performance implications.
In particular, the computation of distance measures and lookback times
are factors of three to four times slower than in the massless neutrino case.
Therefore, if you need to compute many distances in such a cosmology and
performance is critical, it is particularly useful to calculate them on
a grid and use interpolation.

The contribution of photons and neutrinos to the total mass-energy density
can be found as a function of redshift::

  >>> from astropy.cosmology import WMAP7   # WMAP 7-year cosmology
  >>> WMAP7.Ogamma0, WMAP7.Onu0  # Current epoch values  # doctest: +FLOAT_CMP
  (4.985694972799396e-05, 3.442154948307989e-05)
  >>> z = [0, 1.0, 2.0]
  >>> WMAP7.Ogamma(z), WMAP7.Onu(z)  # doctest: +FLOAT_CMP
  (array([4.98603986e-05, 2.74593395e-04, 4.99915942e-04]),
   array([3.44239306e-05, 1.89580995e-04, 3.45145089e-04]))

If you want to exclude photons and neutrinos from your calculations, you can
set ``Tcmb0`` to 0 (which is also the default)::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> import astropy.units as u
  >>> cos = FlatLambdaCDM(70.4 * u.km / u.s / u.Mpc, 0.272, Tcmb0 = 0.0 * u.K)
  >>> cos.Ogamma0, cos.Onu0
  (0.0, 0.0)

You can include photons but exclude any contributions from neutrinos by
setting ``Tcmb0`` to be non-zero (2.725 K is the standard value for our
Universe) but setting ``Neff`` to 0::

  >>> from astropy.cosmology import FlatLambdaCDM
  >>> cos = FlatLambdaCDM(70.4, 0.272, Tcmb0=2.725, Neff=0)
  >>> cos.Ogamma([0, 1, 2])  # Photons are still present  # doctest: +FLOAT_CMP
  array([4.98603986e-05, 2.74642208e-04, 5.00086413e-04])
  >>> cos.Onu([0, 1, 2])  # But not neutrinos  # doctest: +FLOAT_CMP
  array([0., 0., 0.])

The number of neutrino species is assumed to be the floor of ``Neff``,
which in the default case is ``Neff=3``. Therefore, if non-zero neutrino masses
are desired, then three masses should be provided. However, if only one
value is provided, all of the species are assumed to have the same mass.
``Neff`` is assumed to be shared equally between each species.

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
  >>> cosmo.Onu([0, 1.0, 15.0])  # doctest: +FLOAT_CMP
  array([0.00327011, 0.00896845, 0.01257946])
  >>> cosmo.Onu(1) * cosmo.critical_density(1)  # doctest: +FLOAT_CMP
  <Quantity 2.444380380370406e-31 g / cm3>

While these examples used `~astropy.cosmology.FlatLambdaCDM`,
the above examples also apply for all of the other cosmology classes.


For Developers: Using `astropy.cosmology` Inside ``astropy``
============================================================

If you are writing code for the ``astropy`` core or an affiliated package,
it is often useful to assume a default cosmology so that the exact
cosmology does not have to be specified every time a function or method
is called. In this case, it is possible to specify a "default"
cosmology.

You can set the default cosmology to a predefined value by using the
"default_cosmology" option in the ``[cosmology.core]`` section of the
configuration file (see :ref:`astropy_config`). Alternatively, you can
use the ``set`` function of `~astropy.cosmology.default_cosmology` to
set a cosmology for the current Python session. If you have not set a
default cosmology using one of the methods described above, then the
cosmology module will default to using the nine-year WMAP parameters.

It is strongly recommended that you use the default cosmology through
the `~astropy.cosmology.default_cosmology` science state object. An
override option can then be provided using something like the
following::

    def myfunc(..., cosmo=None):
	from astropy.cosmology import default_cosmology

	if cosmo is None:
	    cosmo = default_cosmology.get()

	... your code here ...

This ensures that all code consistently uses the default cosmology
unless explicitly overridden.

.. note::
    In general it is better to use an explicit cosmology (for example
    ``WMAP9.H(0)`` instead of
    ``cosmology.default_cosmology.get().H(0)``). Use of the default
    cosmology should generally be reserved for code that will be
    included in the ``astropy`` core or an affiliated package.

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

See Also
========

* Hogg, "Distance measures in cosmology",
  https://arxiv.org/abs/astro-ph/9905116
* Linder, "Exploring the Expansion History of the Universe", https://arxiv.org/abs/astro-ph/0208512
* NASA's Legacy Archive for Microwave Background Data Analysis,
  https://lambda.gsfc.nasa.gov/

Range of Validity and Reliability
=================================

The code in this sub-package is tested against several widely used
online cosmology calculators and has been used to perform many
calculations in refereed papers. You can check the range of redshifts
over which the code is regularly tested in the module
``astropy.cosmology.tests.test_cosmology``. If you find any bugs,
please let us know by `opening an issue at the GitHub repository
<https://github.com/astropy/astropy/issues>`_!

A more difficult question is the range of redshifts over which
the code is expected to return valid results. This is necessarily
model-dependent, but in general you should not expect the numeric
results to be well behaved for redshifts more than a few times
larger than the epoch of matter-radiation equality (so, for typical
models, not above z = 5-6,000, but for some models much lower redshifts
may be ill-behaved). In particular, you should pay attention to warnings
from the ``scipy`` integration package about integrals failing to converge
(which may only be issued once per session).

The built-in cosmologies use the parameters as listed in the
respective papers. These provide only a limited range of precision,
and so you should not expect derived quantities to match beyond
that precision. For example, the Planck 2013 and 2015 results only provide the
Hubble constant to four digits. Therefore, they should not be expected
to match the age quoted by the Planck team to better than that, despite
the fact that five digits are quoted in the papers.

Reference/API
=============

.. automodapi:: astropy.cosmology
