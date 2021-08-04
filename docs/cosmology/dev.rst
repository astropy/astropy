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

