************
Known Issues
************

.. contents::
   :local:
   :depth: 2

While most bugs and issues are managed using the `astropy issue
tracker <https://github.com/astropy/astropy/issues>`_, this document
lists issues that are too difficult to fix, may require some
intervention from the user to work around, or are caused by bugs in other
projects or packages.

Issues listed on this page are grouped into two categories: The first is known
issues and shortcomings in actual algorithms and interfaces that currently do
not have fixes or workarounds, and that users should be aware of when writing
code that uses ``astropy``. Some of those issues are still platform-specific,
while others are very general. The second category is of common issues that come
up when configuring, building, or installing ``astropy``. This also includes
cases where the test suite can report false negatives depending on the context/
platform on which it was run.

Known Deficiencies
==================

.. _quantity_issues:

Quantities Lose Their Units with Some Operations
------------------------------------------------

Quantities are subclassed from ``numpy``'s `~numpy.ndarray` and while we have
ensured that ``numpy`` functions will work well with them, they do not always
work in functions from ``scipy`` or other packages that use ``numpy``
internally, but ignore the subclass. Furthermore, at a few places in ``numpy``
itself we cannot control the behaviour. For instance, care must be taken when
setting array slices using Quantities::

    >>> import astropy.units as u
    >>> import numpy as np
    >>> a = np.ones(4)
    >>> a[2:3] = 2*u.kg
    >>> a # doctest: +FLOAT_CMP
    array([1., 1., 2., 1.])

::

    >>> a = np.ones(4)
    >>> a[2:3] = 1*u.cm/u.m
    >>> a # doctest: +FLOAT_CMP
    array([1., 1., 1., 1.])

Either set single array entries or use lists of Quantities::

    >>> a = np.ones(4)
    >>> a[2] = 1*u.cm/u.m
    >>> a # doctest: +FLOAT_CMP
    array([1.  , 1.  , 0.01, 1.  ])

::

    >>> a = np.ones(4)
    >>> a[2:3] = [1*u.cm/u.m]
    >>> a # doctest: +FLOAT_CMP
    array([1.  , 1.  , 0.01, 1.  ])

Both will throw an exception if units do not cancel, e.g.::

    >>> a = np.ones(4)
    >>> a[2] = 1*u.cm
    Traceback (most recent call last):
    ...
    TypeError: only dimensionless scalar quantities can be converted to Python scalars


See: https://github.com/astropy/astropy/issues/7582

Multiplying a `pandas.Series` with an `~astropy.units.Unit` does not produce a |Quantity|
-----------------------------------------------------------------------------------------

Quantities may work with certain operations on `~pandas.Series` but
this behaviour is not tested.
For example, multiplying a `~pandas.Series` instance
with a unit will *not* return a |Quantity|. It will return a `~pandas.Series`
object without any unit:

.. doctest-requires:: pandas>=2.0

   >>> import pandas as pd
   >>> import astropy.units as u
   >>> a = pd.Series([1., 2., 3.])
   >>> a * u.m
   0    1.0
   1    2.0
   2    3.0
   dtype: float64

To avoid this, it is best to initialize the |Quantity| directly:

.. doctest-requires:: pandas>=2.0

    >>> u.Quantity(a, u.m)
    <Quantity [1., 2., 3.] m>

Note that the overrides pandas provides are not complete, and
as a consequence, using the (in-place) shift operator does work:

.. doctest-requires:: pandas>=2.0

   >>> b = a << u.m
   >>> b
   <Quantity [1., 2., 3.] m>
   >>> a <<= u.m
   >>> a
   <Quantity [1., 2., 3.] m>

But this is fragile as this may stop working in future versions of
pandas if they decide to override the dunder methods.

See: https://github.com/astropy/astropy/issues/11247

Using Numpy array creation functions to initialize Quantity
-----------------------------------------------------------
Trying the following example will ignore the unit:

    >>> np.full(10, 1 * u.m)
    array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

However, the following works as one would expect

    >>> np.full(10, 1.0, like=u.Quantity([], u.m))
    <Quantity [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.] m>

and is equivalent to::

    >>> np.full(10, 1) << u.m
    <Quantity [1., 1., 1., 1., 1., 1., 1., 1., 1., 1.] m>

`~numpy.zeros`, `~numpy.ones`, and `~numpy.empty` behave similarly.

`~numpy.arange` also supports the ``like`` keyword argument

    >>> np.arange(0 * u.cm, 1 * u.cm, 1 * u.mm, like=u.Quantity([], u.cm))
    <Quantity [0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] cm>

Also note that the unit of the output array is dictated by that of the ``stop``
argument, and that, like for quantities generally, the data has a floating-point
dtype. If ``stop`` is a pure number, the unit of the output will default to that
of the ``like`` argument.

As with ``~numpy.full`` and similar functions, one may alternatively move the
units outside of the call to `~numpy.arange`::

    >>> np.arange(0, 10, 1) << u.mm
    <Quantity [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.] mm>

Or use `~numpy.linspace`:

    >>> np.linspace(0 * u.cm, 9 * u.mm, 10)
    <Quantity [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.] mm>


Quantities Lose Their Units When Broadcasted
--------------------------------------------

When broadcasting Quantities, it is necessary to pass ``subok=True`` to
`~numpy.broadcast_to`, or else a bare `~numpy.ndarray` will be returned::

   >>> q = u.Quantity(np.arange(10.), u.m)
   >>> b = np.broadcast_to(q, (2, len(q)))
   >>> b # doctest: +FLOAT_CMP
   array([[0., 1., 2., 3., 4., 5., 6., 7., 8., 9.],
          [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]])
   >>> b2 = np.broadcast_to(q, (2, len(q)), subok=True)
   >>> b2 # doctest: +FLOAT_CMP
   <Quantity [[0., 1., 2., 3., 4., 5., 6., 7., 8., 9.],
              [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]] m>

This is analogous to the case of passing a Quantity to `~numpy.array`::

   >>> a = np.array(q)
   >>> a # doctest: +FLOAT_CMP
   array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])
   >>> a2 = np.array(q, subok=True)
   >>> a2 # doctest: +FLOAT_CMP
   <Quantity [0., 1., 2., 3., 4., 5., 6., 7., 8., 9.] m>

See: https://github.com/astropy/astropy/issues/7832

Chained Quantity comparisons to dimensionless zero can be misleading
--------------------------------------------------------------------

When chaining comparisons using Quantities and dimensionless zero,
the result may be misleading::

   >>> 0 * u.Celsius == 0 * u.m  # Correct
   False
   >>> 0 * u.Celsius == 0 == 0 * u.m  # Misleading
   np.True_

What the second comparison is really doing is this::

   >>> (0 * u.Celsius == 0) and (0 == 0 * u.m)
   np.True_

See: https://github.com/astropy/astropy/issues/15103

numpy.prod cannot be applied to Quantity
----------------------------------------

Using ``numpy.prod`` function on a Quantity would result in error.
This is because correctly implementing it for Quantity is fairly
difficult, since, unlike for most numpy functions, the result unit
depends on the shape of the input (rather than only on the units
of the inputs).

    >>> np.prod([1, 2, 3] * u.m)
    Traceback (most recent call last):
    ...
    astropy.units.errors.UnitsError: Cannot use 'reduce' method on ufunc multiply with a Quantity instance as it would change the unit.

See: https://github.com/astropy/astropy/issues/18429

def_unit should not be used for logarithmic unit
------------------------------------------------

When defining custom unit involving logarithmic unit, ``def_unit`` usage
should be avoided because it might result in surprising behavior::

    >>> dBW = u.def_unit('dBW', u.dB(u.W))
    >>> 1 * dBW
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand type(s) for *: 'int' and 'Unit'

Instead, it could be defined directly as such::

    >>> dBW = u.dB(u.W)
    >>> 1 * dBW
    <Decibel 1. dB(W)>

See: https://github.com/astropy/astropy/issues/5945

mmap Support for ``astropy.io.fits`` on GNU Hurd
------------------------------------------------

On Hurd and possibly other platforms, ``flush()`` on memory-mapped files are not
implemented, so writing changes to a mmap'd FITS file may not be reliable and is
thus disabled. Attempting to open a FITS file in writeable mode with mmap will
result in a warning (and mmap will be disabled on the file automatically).

See: https://github.com/astropy/astropy/issues/968


Color Printing on Windows
-------------------------

Colored printing of log messages and other colored text does work in Windows,
but only when running in the IPython console. Colors are not currently
supported in the basic Python command-line interpreter on Windows.

``numpy.int64`` does not decompose input ``Quantity`` objects
-------------------------------------------------------------

Python's ``int()`` goes through ``__index__``
while ``numpy.int64`` or ``numpy.int_`` do not go through ``__index__``. This
means that an upstream fix in NumPy is required in order for
``astropy.units`` to control decomposing the input in these functions::

    >>> np.int64((15 * u.km) / (15 * u.imperial.foot))
    np.int64(1)
    >>> np.int_((15 * u.km) / (15 * u.imperial.foot))
    np.int64(1)
    >>> int((15 * u.km) / (15 * u.imperial.foot))
    3280

To convert a dimensionless `~astropy.units.Quantity` to an integer, it is
therefore recommended to use ``int(...)``.

Build/Installation/Test Issues
==============================

Anaconda Users Should Upgrade with ``conda``, Not ``pip``
---------------------------------------------------------

Upgrading ``astropy`` in the Anaconda Python distribution using ``pip`` can result
in a corrupted install with a mix of files from the old version and the new
version. Anaconda users should update with ``conda update astropy``. There
may be a brief delay between the release of ``astropy`` on PyPI and its release
via the ``conda`` package manager; users can check the availability of new
versions with ``conda search astropy``.


Locale Errors in MacOS X and Linux
----------------------------------

On MacOS X, you may see the following error when running ``pip``::

    ...
    ValueError: unknown locale: UTF-8

This is due to the ``LC_CTYPE`` environment variable being incorrectly set to
``UTF-8`` by default, which is not a valid locale setting.

On MacOS X or Linux (or other platforms) you may also encounter the following
error::

    ...
      stderr = stderr.decode(stdio_encoding)
    TypeError: decode() argument 1 must be str, not None

This also indicates that your locale is not set correctly.

To fix either of these issues, set this environment variable, as well as the
``LANG`` and ``LC_ALL`` environment variables to e.g. ``en_US.UTF-8`` using, in
the case of ``bash``::

    export LANG="en_US.UTF-8"
    export LC_ALL="en_US.UTF-8"
    export LC_CTYPE="en_US.UTF-8"

To avoid any issues in future, you should add this line to your e.g.
``~/.bash_profile`` or ``.bashrc`` file.

To test these changes, open a new terminal and type ``locale``, and you should
see something like::

    $ locale
    LANG="en_US.UTF-8"
    LC_COLLATE="en_US.UTF-8"
    LC_CTYPE="en_US.UTF-8"
    LC_MESSAGES="en_US.UTF-8"
    LC_MONETARY="en_US.UTF-8"
    LC_NUMERIC="en_US.UTF-8"
    LC_TIME="en_US.UTF-8"
    LC_ALL="en_US.UTF-8"

If so, you can go ahead and try running ``pip`` again (in the new
terminal).
