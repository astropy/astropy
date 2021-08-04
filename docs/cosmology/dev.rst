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


Astropy Interoperability: I/O and your Cosmology Package
========================================================

If you are developing a package and want to be able to interoperate with
`astropy.cosmology.Cosmology`, you're in the right place! Here we will discuss
and provide examples for enabling Astropy to read and write your file formats,
but also convert your cosmology objects to and from Astropy's |Cosmology|.

The following presumes knowledge of how Astropy structures I/O functions. For
a quick tutorial see :ref:`read_write_cosmologies`.

Now that we know how to build and register functions into
:meth:`~astropy.cosmology.Cosmology.read`,
:meth:`~astropy.cosmology.Cosmology.write`,
:meth:`~astropy.cosmology.Cosmology.from_format`,
:meth:`~astropy.cosmology.Cosmology.to_format`, we can do this in your package.

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

We want to enable conversion between cosmology objects from ``mypckage``
to/from |Cosmology|. All the Astropy interface code is defined in
``mypackage/io/astropy_convert.py``. The following is a rough outline of the
necessary functions and how to register them with astropy's unified
I/O to be automatically available to
:meth:`astropy.cosmology.Cosmology.from_format` and
:meth:`astropy.cosmology.Cosmology.to_format`.

.. literalinclude:: dev/astropy_convert.py
   :language: python


Reading and Writing
-------------------

Everything Astropy read/write related is defined in
``mypackage/io/astropy_io.py``. The following is a rough outline of the read,
write, and identify functions and how to register them with astropy's unified
io to be automatically available to :meth:`astropy.cosmology.Cosmology.read`
and :meth:`astropy.cosmology.Cosmology.write`.

.. literalinclude:: dev/astropy_io.py
   :language: python


If Astropy is an optional dependency
------------------------------------

The ``astropy_io`` and ``astropy_convert`` modules are written assuming Astropy
is installed. If in ``mypackage`` it is an optional dependency then it is
important to detect if Astropy is installed (and the correct version) before
importing ``astropy_io`` and ``astropy_convert``.
We do this in ``mypackage/io/__init__.py``:

.. literalinclude:: dev/io_init.py
   :language: python


Astropy Interoperability Tests
------------------------------

Lastly, it's important to test that everything works. In this example package
all such tests are contained in ``mypackage/io/tests/test_astropy_io.py``.
These tests require Astropy and will be skipped if it is not installed (and
not the correct version), so at least one test in the test matrix should
include ``astropy >= 5.0``.

.. literalinclude:: dev/test_astropy_convert.py
   :language: python

.. literalinclude:: dev/test_astropy_io.py
   :language: python
