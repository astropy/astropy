.. _astropy-cosmology-for-developers:

Cosmology For Developers
************************

Cosmologies in Functions
========================

It is often useful to assume a default cosmology so that the exact cosmology
does not have to be specified every time a function or method is called. In
this case, it is possible to specify a "default" cosmology.

You can set the default cosmology to a predefined value by using the
"default_cosmology" option in the ``[cosmology.core]`` section of the
configuration file (see :ref:`astropy_config`). Alternatively, you can use the
:meth:`~astropy.cosmology.default_cosmology.set` function of
|default_cosmology| to set a cosmology for the current Python session. If you
have not set a default cosmology using one of the methods described above, then
the cosmology module will default to using the
``default_cosmology._value`` parameters.

You can override the default cosmology through the |default_cosmology| science
state object, using something like the following:

.. code-block:: python

    from astropy.cosmology import default_cosmology

    def myfunc(..., cosmo=None):
        if cosmo is None:
            cosmo = default_cosmology.get()

        ... function code here ...

This ensures that all code consistently uses the default cosmology unless
explicitly overridden.

.. note::

    If you are preparing a paper and thus need to ensure your code provides
    reproducible results, it is better to use an explicit cosmology (for
    example ``WMAP9.H(0)`` instead of ``default_cosmology.get().H(0)``).
    Use of the default cosmology should generally be reserved for code that
    allows for the global cosmology state to be changed; e.g. code included in
    ``astropy`` core or an affiliated package.


.. _astropy-cosmology-custom:

Custom Cosmologies
==================

In :mod:`astropy.cosmology` cosmologies are classes, so custom cosmologies may
be implemented by subclassing |Cosmology| (or more likely |FLRW|) and adding
details specific to that cosmology. Here we will review some of those details
and tips and tricks to building a performant cosmology class.

.. code-block:: python

    from astropy.cosmology import FLRW

    @dataclass(frozen=True, repr=False, eq=False)
    class CustomCosmology(FLRW):
        ...  # [details discussed below]


.. _astropy-cosmology-custom-parameters:

Parameters
----------

.. |Parameter| replace:: :class:`~astropy.cosmology.Parameter`

An astropy |Cosmology| is characterized by 1) its class, which encodes the
physics, and 2) its free parameter(s), which specify a cosmological realization.
When defining the former, all parameters must be declared using |Parameter| and
should have values assigned at instantiation.

A |Parameter| is a `descriptor <https://docs.python.org/3/howto/descriptor.html>`_.
When accessed from a class it transparently stores information, like the units
and accepted equivalencies, that might be opaquely contained in the constructor
signature or more deeply in the code. On a cosmology instance, the descriptor
will return the parameter value.

There are a number of best practices. For a reference, this is excerpted from
the definition of |FLRW|.

.. code-block:: python

    @dataclass(frozen=True, repr=False, eq=False)
    class FLRW(Cosmology):

        H0: Parameter = Parameter(doc="Hubble constant as an `~astropy.units.Quantity` at z=0",
                                  unit="km/(s Mpc)", fvalidate="scalar")
        Om0: Parameter = Parameter(doc="Omega matter; matter density/critical density at z=0",
                                   fvalidate="non-negative")
        Ode0: Parameter = Parameter(doc="Omega dark energy; dark energy density/critical density at z=0.",
                                    fvalidate="float")
        Tcmb0: Parameter = Parameter(doc="Temperature of the CMB as `~astropy.units.Quantity` at z=0.",
                                     default=0.0 * u.K, unit="Kelvin", fmt="0.4g", fvalidate="scalar")
        Neff: Parameter = Parameter(doc="Number of effective neutrino species.",
                                    default=3.04, fvalidate="non-negative")
        m_nu: Parameter = Parameter(doc="Mass of neutrino species.",
                                    default=0.0*u.eV, unit="eV", equivalencies=u.mass_energy(), fmt="")
        Ob0: Parameter = Parameter(doc="Omega baryon; baryonic matter density/critical density at z=0.",
                                   default=None)

        @Ob0.validator
        def Ob0(self, param, value):
            """Validate baryon density to None or positive float > matter density."""
            if value is None:
                return value
            value = _validate_non_negative(self, param, value)
            if value > self.Om0:
                raise ValueError("baryonic density can not be larger than total matter density.")
            return value

First note that all the parameters are also arguments in ``__init__()``. This is not
strictly necessary, but is good practice. If the parameter has units (and related
equivalencies) these must be specified on the |Parameter|, as seen in
The "H0" item in :attr:`~astropy.cosmology.FLRW.parameters`.

The next important thing to note is how the parameter value is set, in
``__init__``. |Parameter| allows for a value to be set once (before
auto-locking), so ``self.H0 = H0`` will use this setter and put the value on
"._H0". The advantage of this method over direct assignment to the private
attribute is the use of validators. |Parameter| allows for custom value
validators, using the method-decorator ``validator``, that can check a value's
validity and modify the value, e.g to assign units. If no custom ``validator``
is specified the default is to check if the |Parameter| has defined units and
if so, return the value as a |Quantity| with those units, using all enabled and
the parameter's unit equivalencies.

The last thing to note is pretty formatting for the |Cosmology|. Each
|Parameter| defaults to the `format specification
<https://docs.python.org/3/library/string.html#formatspec>`_ ".3g", but this
may be overridden, like :attr:`~astropy.cosmology.FLRW.Tcmb0` does.

If a new cosmology modifies an existing Parameter, then the
:meth:`~astropy.cosmology.Parameter.clone` method is useful to deep-copy the
parameter and change any constructor argument. For example, see
``FlatFLRWMixin`` in ``astropy.cosmology.flrw`` (also shown below).

.. code-block:: python

    @dataclass(frozen=True, repr=False, eq=False)
    class FlatFLRWMixin(FlatCosmologyMixin):
        ...

        Ode0: Parameter = FLRW.parameters["Ode0"].clone(derived=True)

Mixins
------

`Mixins <https://en.wikipedia.org/wiki/Mixin>`_ are used in
:mod:`~astropy.cosmology` to reuse code across multiple classes in different
inheritance lines. We use the term loosely as mixins are meant to be strictly
orthogonal, but may not be, particularly in ``__init__``.

Currently the only mixin is |FlatCosmologyMixin| and its |FLRW|-specific
subclass |FlatFLRWMixin|. "Flat" cosmologies should use this mixin.
|FlatFLRWMixin| must precede the base class in the multiple-inheritance so that
this mixin's ``__init__`` proceeds the base class'.


.. _astropy-cosmology-fast-integrals:

Speeding up Integrals in Custom Cosmologies
-------------------------------------------

The supplied cosmology classes use a few tricks to speed up distance and time
integrals.  It is not necessary for anyone subclassing |FLRW| to use these
tricks -- but if they do, such calculations may be a lot faster.

The first, more basic, idea is that, in many cases, it's a big deal to provide
explicit formulae for :meth:`~astropy.cosmology.FLRW.inv_efunc` rather than
simply setting up ``de_energy_scale`` -- assuming there is a nice expression.
As noted above, almost all of the provided classes do this, and that template
can pretty much be followed directly with the appropriate formula changes.

The second, and more advanced, option is to also explicitly provide a scalar
only version of :meth:`~astropy.cosmology.FLRW.inv_efunc`. This results in a
fairly large speedup (>10x in most cases) in the distance and age integrals,
even if only done in python, because testing whether the inputs are iterable or
pure scalars turns out to be rather expensive. To take advantage of this, the
key thing is to explicitly set the instance variables
``self._inv_efunc_scalar`` and ``self._inv_efunc_scalar_args`` in the
constructor for the subclass, where the latter are all the arguments except
``z`` to ``_inv_efunc_scalar``. The provided classes do use this optimization,
and in fact go even further and provide optimizations for no radiation, and for
radiation with massless neutrinos coded in cython. Consult the |FLRW|
subclasses and ``scalar_inv_efuncs`` for the details.

However, the important point is that it is *not* necessary to do this.

.. _cosmology_mypackage:

Astropy Interoperability: I/O and your Cosmology Package
========================================================

If you are developing a package and want to be able to interoperate with
|Cosmology|, you're in the right place! Here we will discuss how to enable
Astropy to read and write your file formats, and convert your cosmology objects
to and from Astropy's |Cosmology|.

The following presumes knowledge of how Astropy structures I/O functions. For
a quick tutorial see :ref:`cosmology_io`.

Now that we know how to build and register functions into |Cosmology.read|,
|Cosmology.write|, |Cosmology.from_format|, |Cosmology.to_format|, we can do
this in your package.

Consider a package -- since this is mine, it's cleverly named ``mypackage`` --
with the following file structure: a module for cosmology codes and a module
for defining related input/output functions. In the cosmology module are
defined cosmology classes and a file format -- ``myformat`` -- and everything
should interoperate with astropy. The tests are done with :mod:`pytest` and are
integrated within the code structure.

.. code-block:: text
    :emphasize-lines: 7,8,9,13,14

    mypackage/
        __init__.py
        cosmology/
            __init__.py
            ...
        io/
            __init__.py
            astropy_convert.py
            astropy_io.py
            ...
            tests/
                __init__.py
                test_astropy_convert.py
                test_astropy_io.py
                ...


Converting Objects Between Packages
-----------------------------------

We want to enable conversion between cosmology objects from ``mypackage``
to/from |Cosmology|. All the Astropy interface code is defined in
``mypackage/io/astropy_convert.py``. The following is a rough outline of the
necessary functions and how to register them with astropy's unified I/O to be
automatically available to |Cosmology.from_format| and |Cosmology.to_format|.


Reading and Writing
-------------------

Everything Astropy read/write related is defined in
``mypackage/io/astropy_io.py``. The following is a rough outline of the read,
write, and identify functions and how to register them with astropy's unified
IO to be automatically available to |Cosmology.read| and |Cosmology.write|.


If Astropy is an optional dependency
------------------------------------

The ``astropy_io`` and ``astropy_convert`` modules are written assuming Astropy
is installed. If in ``mypackage`` it is an optional dependency then it is
important to detect if Astropy is installed (and the correct version) before
importing ``astropy_io`` and ``astropy_convert``.
We do this in ``mypackage/io/__init__.py``:


Astropy Interoperability Tests
------------------------------

Lastly, it's important to test that everything works. In this example package
all such tests are contained in ``mypackage/io/tests/test_astropy_io.py``.
These tests require Astropy and will be skipped if it is not installed (and
not the correct version), so at least one test in the test matrix should
include ``astropy >= 5.0``.
