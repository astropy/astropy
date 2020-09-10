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

Quantities are subclassed from NumPy's `~numpy.ndarray` and in some NumPy
operations (and in SciPy operations using NumPy internally) the subclass is
ignored, which means that either a plain array is returned, or a
`~astropy.units.quantity.Quantity` without units.
E.g., prior to astropy 4.0 and numpy 1.17::

    >>> import astropy.units as u
    >>> import numpy as np
    >>> q = u.Quantity(np.arange(10.), u.m)
    >>> np.dot(q,q) # doctest: +SKIP
    285.0
    >>> np.hstack((q,q)) # doctest: +SKIP
    <Quantity [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 0., 1., 2., 3., 4., 5.,
               6., 7., 8., 9.] (Unit not initialised)>

And for all versions::

    >>> ratio = (3600 * u.s) / (1 * u.h)
    >>> ratio # doctest: +FLOAT_CMP
    <Quantity 3600. s / h>
    >>> np.array(ratio) # doctest: +FLOAT_CMP
    array(3600.)
    >>> np.array([ratio]) # doctest: +FLOAT_CMP
    array([1.])

Workarounds are available for some cases. For the above::

    >>> q.dot(q) # doctest: +FLOAT_CMP
    <Quantity 285. m2>

    >>> np.array(ratio.to(u.dimensionless_unscaled)) # doctest: +FLOAT_CMP
    array(1.)

    >>> u.Quantity([q, q]).flatten() # doctest: +FLOAT_CMP
    <Quantity [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 0., 1., 2., 3., 4., 5.,
               6., 7., 8., 9.] m>

An incomplete list of specific functions which are known to exhibit
this behavior (prior to astropy 4.0 and numpy 1.17) follows:

* `numpy.dot`
* `numpy.hstack`, `numpy.vstack`, ``numpy.c_``, ``numpy.r_``, `numpy.append`
* `numpy.where`
* `numpy.choose`
* `numpy.vectorize`
* pandas DataFrame(s)


See: https://github.com/astropy/astropy/issues/1274


Care must be taken when setting array slices using Quantities::

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
    >>> a[2] = 1*u.cm # doctest: +SKIP
    Traceback (most recent call last):
    ...
    TypeError: only dimensionless scalar quantities can be converted to Python scalars


See: https://github.com/astropy/astropy/issues/7582

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

Quantities Float Comparison with np.isclose Fails
-------------------------------------------------

Comparing Quantities floats using the NumPy function `~numpy.isclose` fails on
NumPy versions before 1.17 as the comparison between ``a`` and ``b``
is made using the formula

.. math::

    |a - b| \le (a_\textrm{tol} + r_\textrm{tol} \times |b|)

This will result in the following traceback when using this with Quantities::

    >>> from astropy import units as u, constants as const
    >>> import numpy as np
    >>> np.isclose(500 * u.km/u.s, 300 * u.km / u.s)  # doctest: +SKIP
    Traceback (most recent call last):
    ...
    UnitConversionError: Can only apply 'add' function to dimensionless quantities when other argument is not a quantity (unless the latter is all zero/infinity/nan)

If one cannot upgrade to numpy 1.17 or later, one solution is::

    >>> np.isclose(500 * u.km/u.s, 300 * u.km / u.s, atol=1e-8 * u.mm / u.s)
    False

Quantities in np.linspace Failure on NumPy 1.10
-----------------------------------------------

`~numpy.linspace` does not work correctly with quantities when using NumPy
1.10.0 to 1.10.5 due to a bug in NumPy. The solution is to upgrade to NumPy
1.10.6 or later, in which the bug was fixed.


mmap Support for ``astropy.io.fits`` on GNU Hurd
------------------------------------------------

On Hurd and possibly other platforms, ``flush()`` on memory-mapped files are not
implemented, so writing changes to a mmap'd FITS file may not be reliable and is
thus disabled. Attempting to open a FITS file in writeable mode with mmap will
result in a warning (and mmap will be disabled on the file automatically).

See: https://github.com/astropy/astropy/issues/968


Bug with Unicode Endianness in ``io.fits`` for Big Endian Processors
--------------------------------------------------------------------

On big endian processors (e.g. SPARC, PowerPC, MIPS), string columns in FITS
files may not be correctly read when using the ``Table.read`` interface. This
will be fixed in a subsequent bug fix release of ``astropy`` (see `bug report here
<https://github.com/astropy/astropy/issues/3415>`_).


Color Printing on Windows
-------------------------

Colored printing of log messages and other colored text does work in Windows,
but only when running in the IPython console. Colors are not currently
supported in the basic Python command-line interpreter on Windows.

``numpy.int64`` does not decompose input ``Quantity`` objects
-------------------------------------------------------------

Python's ``int()`` goes through ``__index__``
while ``numpy.int64`` or ``numpy.int_`` do not go through ``__index__``. This
means that an upstream fix in ``numpy` is required in order for
``astropy.units`` to control decomposing the input in these functions::

    >>> np.int64((15 * u.km) / (15 * u.imperial.foot))
    1
    >>> np.int_((15 * u.km) / (15 * u.imperial.foot))
    1
    >>> int((15 * u.km) / (15 * u.imperial.foot))
    3280

To convert a dimensionless `~astropy.units.Quantity` to an integer, it is
therefore recommended to use ``int(...)``.

Inconsistent behavior when converting complex numbers to floats
---------------------------------------------------------------

Attempting to use `float` or NumPy's ``numpy.float`` on a standard
complex number (e.g., ``5 + 6j``) results in a `TypeError`.  In
contrast, using `float` or ``numpy.float`` on a complex number from
NumPy (e.g., ``numpy.complex128``) drops the imaginary component and
issues a ``numpy.ComplexWarning``.  This inconsistency persists between
`~astropy.units.Quantity` instances based on standard and NumPy
complex numbers.  To get the real part of a complex number, it is
recommended to use ``numpy.real``.

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

On MacOS X, you may see the following error when running ``setup.py``::

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

If so, you can go ahead and try running ``setup.py`` again (in the new
terminal).


Creating a Time Object Fails with ValueError After Upgrading ``astropy``
------------------------------------------------------------------------

In some cases, when users have upgraded ``astropy`` from an older version to v1.0
or greater, they have run into the following crash when trying to create an
`~astropy.time.Time` object::

    >>> from astropy.time import Time
    >>> datetime = Time('2012-03-01T13:08:00', scale='utc') # doctest: +SKIP
    Traceback (most recent call last):
    ...
    ValueError: Input values did not match any of the formats where
    the format keyword is optional [u'astropy_time', u'datetime',
    u'jyear_str', u'iso', u'isot', u'yday', u'byear_str']

This problem can occur when there is a version mismatch between the compiled
ERFA library (included as part of ``astropy`` in most distributions), and
the version of the ``astropy`` Python source.

This can be from a number of causes. The most likely is that when installing the
new ``astropy`` version, your previous ``astropy`` version was not fully uninstalled
first, resulting in a mishmash of versions. Your best bet is to fully remove
``astropy`` from its installation path and reinstall from scratch using your
preferred installation method. Removing the old version may be achieved by
removing the entire ``astropy/`` directory from within the
``site-packages`` directory it is installed in. However, if in doubt, ask
how best to uninstall packages from your preferred Python distribution.

Another possible cause of this error, in particular for people developing on
Astropy and installing from a source checkout, is that your Astropy build
directory is unclean. To fix this, run ``git clean -dfx``. This removes
*all* build artifacts from the repository that aren't normally tracked by git.
Make sure before running this that there are no untracked files in the
repository you intend to save. Then rebuild/reinstall from the clean repo.


Failing Logging Tests When Running the Tests in IPython
-------------------------------------------------------

When running the Astropy tests using ``astropy.test()`` in an IPython
interpreter, some of the tests in the ``astropy/tests/test_logger.py`` *might*
fail depending on the version of IPython or other factors.
This is due to mutually incompatible behaviors in IPython and pytest, and is
not due to a problem with the test itself or the feature being tested.

See: https://github.com/astropy/astropy/issues/717


Some Docstrings Can Not Be Displayed in IPython < 0.13.2
--------------------------------------------------------

Displaying long docstrings that contain Unicode characters may fail on
some platforms in the IPython console (prior to IPython version
0.13.2)::

    In [1]: import astropy.units as u

    In [2]: u.Angstrom?
    Out[2]: ERROR: UnicodeEncodeError: 'ascii' codec can't encode character u'\xe5' in
    position 184: ordinal not in range(128) [IPython.core.page]

This can be worked around by changing the default encoding to ``utf-8``
by adding the following to your ``sitecustomize.py`` file::

    import sys
    sys.setdefaultencoding('utf-8')

Note that in general, `this is not recommended
<https://stackoverflow.com/questions/3828723/why-should-we-not-use-sys-setdefaultencodingutf-8-in-a-py-script>`_,
because it can hide other Unicode encoding bugs in your application.
However, if your application does not deal with text
processing and you just want docstrings to work, this may be
acceptable.

The IPython issue: https://github.com/ipython/ipython/pull/2738

Compatibility Issues with pytest 3.7 and later
----------------------------------------------

Due to a bug in `pytest <http://www.pytest.org>`_ related to test collection,
the tests for the core ``astropy`` package for version 2.0.x (LTS), and for
packages using the core package's test infrastructure and being tested against
2.0.x (LTS), will not be executed correctly with pytest 3.7, 3.8, or 3.9. The
symptom of this bug is that no tests or only tests in RST files are collected.
In addition, ``astropy`` 2.0.x (LTS) is not compatible with pytest 4.0 and above,
as in this case deprecation errors from pytest can cause tests to fail.
Therefore, when testing against ``astropy`` v2.0.x (LTS), pytest 3.6 or earlier
versions should be used. These issues do not occur in version 3.0.x and above of
the core package.

There is an unrelated issue that also affects more recent versions of
``astropy`` when testing with pytest 4.0 and later, which can
cause issues when collecting tests — in this case, the symptom is that the
test collection hangs and/or appears to run the tests recursively. If you are
maintaining a package that was created using the Astropy
`package template <https://github.com/astropy/package-template>`_, then
this can be fixed by updating to the latest version of the ``_astropy_init.py``
file. The root cause of this issue is that pytest now tries to pick up the
top-level ``test()`` function as a test, so we need to make sure that we set a
``test.__test__`` attribute on the function to ``False``.
