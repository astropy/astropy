.. _astropy-cosmology-for-developers:

Cosmology For Developers
************************

Cosmologies in Functions
========================

It is often useful to assume a default cosmology so that the exact
cosmology does not have to be specified every time a function or method
is called. In this case, it is possible to specify a "default"
cosmology.

You can set the default cosmology to a predefined value by using the
"default_cosmology" option in the ``[cosmology.core]`` section of the
configuration file (see :ref:`astropy_config`). Alternatively, you can
use the :meth:`~astropy.cosmology.default_cosmology.set` function of
:func:`~astropy.cosmology.default_cosmology` to set a cosmology for the current
Python session. If you have not set a default cosmology using one of the
methods described above, then the cosmology module will default to using the
``default_cosmology._value`` parameters.

It is strongly recommended that you use the default cosmology through
the :class:`~astropy.cosmology.default_cosmology` science state object. An
override option can then be provided using something like the
following:

.. code-block:: python

    from astropy.cosmology import default_cosmology

    def myfunc(..., cosmo=None):
        if cosmo is None:
            cosmo = default_cosmology.get()

        ... function code here ...

This ensures that all code consistently uses the default cosmology
unless explicitly overridden.

.. note::

    In general it is better to use an explicit cosmology (for example
    ``WMAP9.H(0)`` instead of
    ``cosmology.default_cosmology.get().H(0)``). Use of the default
    cosmology should generally be reserved for code that will be
    included in the ``astropy`` core or an affiliated package.


.. _astropy-cosmology-fast-integrals:

Speeding up Integrals in Custom Cosmologies
===========================================

The supplied cosmology classes use a few tricks to speed up distance and time
integrals.  It is not necessary for anyone subclassing
:class:`~astropy.cosmology.FLRW` to use these tricks -- but if they do, such
calculations may be a lot faster.

The first, more basic, idea is that, in many
cases, it's a big deal to provide explicit formulae for
:meth:`~astropy.cosmology.FLRW.inv_efunc` rather than simply setting up
``de_energy_scale`` -- assuming there is a nice expression. As noted above,
almost all of the provided classes do this, and that template can pretty much
be followed directly with the appropriate formula changes.

The second, and more advanced, option is to also explicitly provide a
scalar only version of :meth:`~astropy.cosmology.FLRW.inv_efunc`. This results
in a fairly large speedup (>10x in most cases) in the distance and age
integrals, even if only done in python, because testing whether the inputs are
iterable or pure scalars turns out to be rather expensive. To take advantage of
this, the key thing is to explicitly set the instance variables
``self._inv_efunc_scalar`` and ``self._inv_efunc_scalar_args`` in the
constructor for the subclass, where the latter are all the arguments except
``z`` to ``_inv_efunc_scalar``. The provided classes do use this optimization,
and in fact go even further and provide optimizations for no radiation, and for
radiation with massless neutrinos coded in cython. Consult the 
:class:`~astropy.cosmology.FLRW` subclasses for details, and
``scalar_inv_efuncs`` for the details.

However, the important point is that it is *not* necessary to do this.

